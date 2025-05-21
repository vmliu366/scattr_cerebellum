# Directories
zona_dir = str(Path(config["output_dir"]) / "zona_bb_subcortex")
labelmerge_dir = str(Path(config["output_dir"]) / "labelmerge")
log_dir = str(Path(config["output_dir"]) / ".logs" / "zona_bb_subcortex")

# Make directory if it doesn't exist
Path(zona_dir).mkdir(parents=True, exist_ok=True)


# BIDS partials
bids_anat = partial(
    bids,
    root=zona_dir,
    datatype="anat",
    **inputs_t1w.subj_wildcards,
)

bids_labelmerge = partial(
    bids,
    root=str(Path(labelmerge_dir) / "combined")
    if not config.get("skip_labelmerge")
    else config.get("labelmerge_base_dir") or zona_dir,
    **inputs_t1w.subj_wildcards,
)

bids_log = partial(
    bids,
    root=log_dir,
    **inputs_t1w.subj_wildcards,
)

""" References:
J.C. Lau, Y. Xiao, R.A.M. Haast, G. Gilmore, K. Uludağ, K.W. MacDougall, 
R.S. Menon, A.G. Parrent, T.M. Peters, A.R. Khan. Direct visualization and 
characterization of the human zona incerta and surrounding structures. 
Hum. Brain Mapp., 41 (2020), pp. 4500-4517, 10.1002/hbm.25137

Y. Xiao, J.C. Lau, T. Anderson, J. DeKraker, D.L. Collins, T. Peters, 
A.R. Khan. An accurate registration of the BigBrain dataset with the MNI PD25 
and ICBM152 atlases. Sci. Data, 6 (2019), p. 210, 10.1038/s41597-019-0217-0
"""


rule cp_zona_tsv:
    """Copy tsv to zona dir"""
    input:
        zona_tsv=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"]["tsv"])
        ),
    output:
        zona_tsv=f"{zona_dir}/desc-ZonaBB_dseg.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        time=10,
    group:
        "dseg_tsv"
    shell:
        "cp -v {input.zona_tsv} {output.zona_tsv}"


rule reg2native:
    """
    Create transforms from chosen template space to subject native space via 
    ANTsRegistrationSyNQuick
    """
    input:
        template=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["T1w"])
        ),
        target=lambda wildcards: inputs_t1w["T1w"]
        .filter(**wildcards)
        .expand()[0],
    params:
        out_dir=directory(str(Path(bids_anat()).parent)),
        out_prefix=bids_anat(
            desc=f"from{config['Space']}toNative_",
        ),
    output:
        t1w_nativespace=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="Warped.nii.gz",
        ),
        warp=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="1Warp.nii.gz",
        ),
        affine=bids_anat(
            desc=f"from{config['Space']}toNative",
            suffix="0GenericAffine.mat",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
        threads=4,
    log:
        bids_log(suffix="reg2native.log"),
    group:
        "subcortical_1"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        echo {input.target}        
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}

        echo "Computed output directory: {params.out_dir}"
        mkdir -p {params.out_dir}

        antsRegistrationSyNQuick.sh -n {threads} -d 3 \\
        -f {input.target} -m {input.template} \\
        -o {params.out_prefix} &> {log}
        """


rule warp2native:
    """Warp subcortical parcellations to subject native space"""
    input:
        dseg=str(
            Path(workflow.basedir).parent
            / Path(config["zona_bb_subcortex"][config["Space"]]["dir"])
            / Path(config["zona_bb_subcortex"][config["Space"]]["seg"])
        ),
        target=rules.reg2native.input.target,
        warp=rules.reg2native.output.warp,
        affine=rules.reg2native.output.affine,
    output:
        nii=bids_anat(
            space="T1w",
            desc="ZonaBB",
            suffix="dseg.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=30,
        threads=4,
    log:
        bids_log(suffix="warp2native.log"),
    group:
        "subcortical_1"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} 

        antsApplyTransforms -v -d 3 -n MultiLabel \\
        -i {input.dseg} -r {input.target} \\
        -t {input.warp} -t {input.affine} \\
        -o {output.nii} &> {log}
        """

# ====================================== cerebellum-specific ========================================

""" References:
Diedrichsen, J., Balsters, J. H., Flavell, J., Cussans, E., & Ramnani, N. (2009). 
A probabilistic atlas of the human cerebellum. Neuroimage.

Diedrichsen, J., Maderwald, S., Kuper, M., Thurling, M., Rabe, K., Gizewski, E. R., et al. (2011). 
Imaging the deep cerebellar nuclei: A probabilistic atlas and normalization procedure. Neuroimage.

http://www.diedrichsenlab.org/imaging/propatlas.htm
https://github.com/DiedrichsenLab/cerebellar_atlases

"""

rule reg2native_cerebellum:
    """
    Create a transform from the Diedrichsen atlas space into subject-native space.
    """
    input:
        template=lambda wc: str(
            Path(workflow.basedir).parent
            / config["cerebellum"][config["Space"]]["T1w"]
        ),
        target=lambda wc: inputs_t1w["T1w"].filter(**wc).expand()[0],
    params:
        out_dir  = directory(str(Path(bids_anat()).parent)),
        out_prefix = bids_anat(desc="fromCerebellumToNative_"),
    output:
        warped_t1w = bids_anat(desc="fromCerebellumToNative", suffix="Warped.nii.gz"),
        warp_mat   = bids_anat(desc="fromCerebellumToNative", suffix="1Warp.nii.gz"),
        affine_mat = bids_anat(desc="fromCerebellumToNative", suffix="0GenericAffine.mat"),
    threads:    4
    resources:
        mem_mb=16000, time=60, threads=4,
    log:
        bids_log(suffix="reg2native_cerebellum.log")
    group:
        "subcortical_1"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        antsRegistrationSyNQuick.sh -n {threads} -d 3 \
            -f {input.target} -m {input.template} \
            -o {params.out_prefix} &> {log}
        """

rule warp2native_cerebellum:
    """
    Warp the cerebellum atlas segmentation into subject-native space.
    """
    input:
        seg = lambda wc: str(Path(workflow.basedir).parent
            / config["cerebellum"][config["Space"]]["seg"]),
        target = rules.reg2native_cerebellum.input.target,
        warp   = rules.reg2native_cerebellum.output.warp_mat,
        affine = rules.reg2native_cerebellum.output.affine_mat,
    output:
        nii = bids_anat(space="T1w", desc="Cerebellum", suffix="dseg.nii.gz"),
    threads: 4
    resources:
        mem_mb=16000, time=30, threads=4,
    log:
        bids_log(suffix="warp2native_cerebellum.log")
    group:
        "subcortical_1"
    container:
        config["singularity"]["scattr"]
    shell:
        """
        antsApplyTransforms \
           -i {input.seg} \
           -r {input.target} \
           -t {input.warp} -t {input.affine} \
           -o {output.nii} &> {log}
        """

rule cp_cerebellum_tsv:
    """Copy cerebellum atlas TSV into the zona_bb_subcortex BIDS output dir"""
    input:
        # repo‐relative path to your source TSV as defined in config.yaml
        cereb_tsv = str(
            Path(workflow.basedir).parent
            / Path(config["cerebellum"]["tsv"])
        ),
    output:
        # drop into the same zona_dir, with desc‑Cerebellum label
        cereb_tsv = f"{zona_dir}/desc-Cerebellum_dseg.tsv",
    threads: 1
    resources:
        mem_mb = 4000,
        time   = 10,
        threads=4,
    group:
        "dseg_tsv"
    shell:
        "cp -v {input.cereb_tsv} {output.cereb_tsv}"


# =============================================================================================

rule labelmerge:
    input:
        zona_seg=inputs_t1w["T1w"].expand(
            rules.warp2native.output.nii, allow_missing=False
        )
        if not config.get("labelmerge_base_dir")
        else [],
        fs_seg=inputs_t1w["T1w"].expand(
            rules.fs_xfm_to_native.output.thal, allow_missing=False
        )
        if not config.get("labelmerge_overlay_dir")
        else [],
        fs_tsv=rules.cp_fs_tsv.output.fs_tsv
        if not config.get("labelmerge_overlay_dir")
        else [],
        zona_tsv=rules.cp_zona_tsv.output.zona_tsv
        if not config.get("labelmerge_base_dir")
        else [],
    params:
        labelmerge_out_dir=directory(labelmerge_dir),
        labelmerge_base_dir=(
            config.get("labelmerge_base_dir")
            if config.get("labelmerge_base_dir")
            else zona_dir
        ),
        base_drops=f"--base_drops {config['labelmerge_base_drops']}",
        base_desc=f"--base_desc {config['labelmerge_base_desc']}",
        base_exceptions=(
            f"--base_exceptions {config.get('labelmerge_base_exceptions')}"
            if config.get("labelmerge_base_exceptions")
            else ""
        ),
        overlay_dir=(
            f"--overlay_bids_dir {config.get('labelmerge_overlay_dir') if config.get('labelmerge_overlay_dir') else rules.thalamic_segmentation.input.freesurfer_dir}"
        ),
        overlay_desc=f"--overlay_desc {config['labelmerge_overlay_desc']}",
        overlay_drops=(
            f"--overlay_drops {config.get('labelmerge_overlay_drops')}"
            if config.get("labelmerge_overlay_drops")
            else ""
        ),
        overlay_exceptions=(
            f"--overlay_exceptions {config.get('labelmerge_overlay_exceptions')}"
            if config.get("labelmerge_overlay_exceptions")
            else ""
        ),
        output_desc="combined",
    output:
        seg=inputs_t1w["T1w"].expand(
            bids_labelmerge(
                space="T1w",
                datatype="anat",
                desc="combined",
                suffix="dseg.nii.gz",
            ),
            allow_missing=True,
        ),
        tsv=inputs_t1w["T1w"].expand(
            bids_labelmerge(
                space="T1w",
                datatype="anat",
                desc="combined",
                suffix="dseg.tsv",
            ),
            allow_missing=True,
        ),
    resources:
        mem_mb=16000,
        time=60,
    group:
        "subcortical_group"
    # container:
    #     config["singularity"]["scattr"]
    shell:
        """
        labelmerge {params.labelmerge_base_dir} {params.labelmerge_out_dir} \\
            participant \\
            {params.base_desc} {params.base_drops} {params.base_exceptions} \\
            {params.overlay_dir} {params.overlay_desc} \\
            {params.overlay_drops} {params.overlay_exceptions} \\
            --output_desc {params.output_desc} \\
            --cores {threads} --forceall
        """

# =============================== chaining labelmerge ==================================
rule merge_with_cerebellum:
    """
    Merge the zonaBB+FS ‘combined’ atlas (desc‑combined) with the cerebellum atlas
    via the stand‑alone labelmerge.py script, to produce desc‑combined2.
    """
    input:
        base_seg=inputs_t1w["T1w"].expand(
            rules.labelmerge.output.seg, allow_missing=False
        ),
        cereb_seg=inputs_t1w["T1w"].expand(
            rules.warp2native_cerebellum.output.nii, allow_missing=False
        ),
        base_tsv=inputs_t1w["T1w"].expand(
            bids_labelmerge(
                space="T1w",
                datatype='anat',
                desc="combined",
                suffix="dseg.tsv",
            ),
        ),
        cereb_tsv=rules.cp_cerebellum_tsv.output.cereb_tsv,
    params:
        base_dir=directory(labelmerge_dir),
        cereb_out_dir=directory(labelmerge_dir),
        base_desc="combined",
        overlay_dir=directory(zona_dir),
        overlay_desc="Cerebellum",
        output_desc="combined2",
    output:
        seg=inputs_t1w["T1w"].expand(
            bids_labelmerge(
                space="T1w",
                datatype="anat",
                desc="combined2",
                suffix="dseg.nii.gz",
            ),
            allow_missing=False,
        ),
        tsv=inputs_t1w["T1w"].expand(
            bids_labelmerge(
                space="T1w",
                datatype="anat",
                desc="combined2",
                suffix="dseg.tsv",
            ),
            allow_missing=False,
        ),
    resources:
        mem_mb=16000,
        time=60,
    group:
        "subcortical_group"
    # container:
    #     config["singularity"]["scattr"]
    shell:
        """
        labelmerge {params.base_dir} {params.cereb_out_dir} \\
            participant \\
            --base-desc {params.base_desc} \\
            --overlay-bids-dir {params.overlay_dir} --overlay_desc {params.overlay_desc} \\
            --output_desc {params.output_desc} \\
            --cores {threads} --forceall
        """

# rule merge_with_cerebellum:
#     """
#     Merge the zonaBB+FS ‘combined’ atlas (desc‑combined) with the cerebellum atlas
#     via the stand‑alone labelmerge.py script, to produce desc‑combined2.

#     Currently calling labelmerge.py directly. TO DO: integrate into labelmerge 
#     snakemake workflow 
#     """
#     input:

#         #base_map=rules.labelmerge.output.seg,
#         #base_tsv=rules.labelmerge.output.tsv,
#         #overlay_map = rules.warp2native_cerebellum.output.nii,
#         #overlay_tsv = rules.cp_cerebellum_tsv.output.cereb_tsv,
        
#         # base_map=lambda wildcards: f"{config['output_dir']}/labelmerge/combined/sub-{wildcards.subject}/sub-{wildcards.subject}_space-T1w_desc-combined_dseg.nii.gz",
#         # base_tsv=lambda wildcards: f"{config['output_dir']}/labelmerge/combined/sub-{wildcards.subject}/sub-{wildcards.subject}_space-T1w_desc-combined_dseg.tsv",
#         # overlay_map=lambda wildcards: f"{config['output_dir']}/zona_bb_subcortex/sub-{wildcards.subject}/anat/sub-{wildcards.subject}_space-T1w_desc-Cerebellum_dseg.nii.gz",
#         # overlay_tsv=lambda wildcards: f"{config['output_dir']}/zona_bb_subcortex/desc-Cerebellum_dseg.tsv",
    
#     output:
#         seg = bids_labelmerge(
#             space="T1w", desc="combined2", suffix="dseg.nii.gz"
#         ),
#         tsv = bids_labelmerge(
#             space="T1w", desc="combined2", suffix="dseg.tsv"
#         ),
#     threads: 1
#     resources:
#         mem_mb=4000,
#         time=10,
#     shell:
#         """
#         # call the pure‐Python merger directly  
#         python3 {workflow.basedir}/scripts/zona_bb_subcortex/labelmerge.py \
#             {input.base_map}  {input.base_tsv} \
#             {input.overlay_map} {input.overlay_tsv} \
#             {output.seg}      {output.tsv}
#         """

# ===================================================================================

rule get_num_nodes:
    input:
       seg=bids_labelmerge(
           space="T1w",
        #    datatype="anat" if config.get("skip_labelmerge") else "",
           datatype="anat",
           desc="combined2"
           if not config.get("skip_labelmerge")
           else config.get("labelmerge_base_desc"),
           suffix="dseg.nii.gz",
       ),
        # seg = rules.merge_with_cerebellum.output.seg,
    output:
        num_labels=bids_labelmerge(
                space="T1w",
                datatype="anat" if config.get("skip_labelmerge") else "",
#                desc="combined"
                desc="combined2"
                if not config.get("skip_labelmerge")
                else config.get("labelmerge_base_desc"),
                suffix="numNodes.txt",
            ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
        threads=4,
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/get_num_labels.py"


rule binarize:
    input:
        seg=rules.get_num_nodes.input.seg,
    output:
        mask=bids_labelmerge(
            space="T1w",
            desc="combined2"
            if not config.get("skip_labelmerge")
            else config.get("labelmerge_base_desc"),
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
        threads=4,
    log:
        bids_log(suffix="binarize.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/create_seg_mask.py"


rule add_brainstem:
    input:
        mask=rules.binarize.output.mask,
        aparcaseg=rules.fs_xfm_to_native.output.aparcaseg
        if not config.get("skip_brainstem")
        else [],
    output:
        mask=bids_labelmerge(
            space="T1w",
            desc="labelmergeStem",
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=10,
        threads=4,
    log:
        bids_log(suffix="addBrainstem.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/add_brainstem.py"


rule create_convex_hull:
    input:
        bin_seg=rules.binarize.output.mask
        if config["skip_brainstem"]
        else rules.add_brainstem.output.mask,
    output:
        convex_hull=bids_labelmerge(
            space="T1w",
            desc="ConvexHull",
            suffix="mask.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
        threads=4,
    log:
        bids_log(suffix="createConvexHull.log"),
    group:
        "subcortical_2"
    container:
        config["singularity"]["scattr"]
    script:
        "../scripts/zona_bb_subcortex/convexHull_roi.py"
