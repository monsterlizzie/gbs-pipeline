process GENERATE_SAMPLE_REPORT {
    label 'python_container'
    label 'farm_low'

    tag "$sample_id"

    publishDir "${params.output}/sample_reports", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}_process_report_?.csv")

    output:
    path "${sample_id}_report.csv", emit: report

    script:
    """
    generate_sample_report.py "$sample_id" "${sample_id}_report.csv" ${sample_id}_process_report_*.csv
    """
}

process GENERATE_OVERALL_REPORT {
    label 'python_container'
    label 'farm_low'

    publishDir "${params.output}", mode: "copy"

    input:
    path '*_report.csv'

    output: 
    path "results.csv", emit: report // for combining typer reports
    //path "$overall_report", emit: report // for combining typer reports

    script:
    //input_pattern='*_report.csv' // for combining typer reports
    //overall_report='results.csv' // for combining typer reports
    """
    generate_overall_report.py *_report.csv results.csv
    """
}
