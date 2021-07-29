configfile: "/storage/scratch/covid/container/config/config.yaml"


home_dir = config["container"]
run = config["current_run"]
scripts_dir = os.path.join(home_dir, "scripts/")

rule all:
    input:
        expand(home_dir+ "pangolin_analysis/{run}/lineage_report_{run}_v3_1.csv",run=run),
        expand(home_dir+ "pangolin_analysis/pangolin_lineage_assignment_for_pipeline_{run}.xlsx",run=run),
        expand(home_dir+ "variant_analysis/final_result_with_var_assignment_{run}.xlsx",run=run),
        expand(home_dir+ "output/{run}/database/4_Curated_MCOV_MRN_Strains.xlsx",run=run),
        expand(home_dir+ "output/{run}/database/5_strains_run_group_analysis_table.csv",run=run),
        expand(home_dir+ "output/{run}/database/5_strains_run_group_analysis_table_with_samplesheet.csv",run=run),
        expand(home_dir+ "output/{run}/database/run_samplesheet_match_.xlsx",run=run),
        expand(home_dir+ "output/{run}/1_genomic_mutations_{run}.xlsx",run=run),
        expand(home_dir+ "output/{run}/1_genomic_mutations_all.xlsx",run=run),
        expand(home_dir+ "output/{run}/3_variant_screening_per_variant_summary_pangolin_{run}.xlsx",run=run),
        expand(home_dir+ "output/{run}/4_mcov_strain_variant_map_covid_pangolin_db_input_{run}.csv",run=run),
        expand(home_dir+ "output/{run}/spike_protein_snp_table_{run}_.xlsx",run=run),
        expand(home_dir+ "output/{run}/variant_growth.csv",run=run),
        expand(home_dir+ "output/{run}/{run}_source_file.xlsx",run=run)


rule pangolin_analysis:
    input:
        script = scripts_dir + "pangolin_analysis.sh"
 
    output:
        out_run = home_dir+ "pangolin_analysis/{run}/lineage_report_{run}_v3_1.csv",
        out_combine = home_dir+ "pangolin_analysis/pangolin_lineage_assignment_for_pipeline_{run}.xlsx"
    
    shell: 
        "bash {input.script} -r {run} "

rule variant_analysis:
    input:
        script = scripts_dir + "variant_analysis_combine.py"

    output:
        out_combine = home_dir + "variant_analysis/final_result_with_var_assignment_{run}.xlsx",
    
    shell: 
        "python {input.script} "

rule source_preprocessing:
    input:
        script = scripts_dir + "run_db_update.py"

    output:
        out_combine = home_dir+ "output/{run}/database/4_Curated_MCOV_MRN_Strains.xlsx",
        out_table = home_dir+ "output/{run}/database/5_strains_run_group_analysis_table.csv",
        out_tablewithss = home_dir+ "output/{run}/database/5_strains_run_group_analysis_table_with_samplesheet.csv",
        out_match = home_dir+ "output/{run}/database/run_samplesheet_match_.xlsx"

    shell: 
        "python {input.script} "

rule genomics:
    input:
        script = scripts_dir + "run_genomics.py",
        processed = rules.source_preprocessing.output.out_tablewithss

    output:
        out_mutations_run = home_dir+ "output/{run}/1_genomic_mutations_{run}.xlsx",

    shell: 
        "python {input.script} "

rule snp_analysis:
    input:
        script = scripts_dir + "run_snp_analysis.py",
        processed_run = rules.genomics.output.out_mutations_run,
        pangolin = rules.pangolin_analysis.output.out_combine,

    output:
        out_mutations = home_dir+ "output/{run}/1_genomic_mutations_all.xlsx",
        out_screening = home_dir+ "output/{run}/3_variant_screening_per_variant_summary_pangolin_{run}.xlsx",
        out_coviddb = home_dir+ "output/{run}/4_mcov_strain_variant_map_covid_pangolin_db_input_{run}.csv",
        out_spike_snvs = home_dir+ "output/{run}/spike_protein_snp_table_{run}_.xlsx"

    shell: 
        "python {input.script} "

rule variant_growth:
    input:
        script = scripts_dir + "variant_growth.py",
        processed = rules.snp_analysis.output.out_coviddb

    params:
        start_date = "1-1-2021",
        end_date = "7-20-2021"

    output:
        out_variant_growth = home_dir+ "output/{run}/variant_growth.csv"

    shell: 
        "python {input.script} {params.start_date}  {params.end_date}"

rule post_process:
    input:
        script = scripts_dir + "post_process.py",
        genomics = rules.snp_analysis.output.out_mutations,
        coviddb = rules.snp_analysis.output.out_coviddb,
        variant = rules.variant_analysis.output.out_combine

    output:
        out_source = home_dir+ "output/{run}/{run}_source_file.xlsx"

    shell: 
        "python {input.script} "

