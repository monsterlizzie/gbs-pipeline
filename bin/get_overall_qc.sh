# Determine overall QC result based on File Validity, Read QC, Assembly QC, Mapping QC and Taxonomy QC
# In case File Validity is not PASS, save its value (i.e. description of the issue) to Overall QC
# In case of preprocess/assembly/mapping/taxonomy failure, there will be no relevant QC input, save corresponding MODULE FAILURE to Overall QC

assign_overall_qc() {
    if [[ "$FILE_VALIDITY" == "null" ]]; then
        OVERALL_QC="FILE VALIDATION FAILURE"
        return
    fi
    
    if [[ ! "$FILE_VALIDITY" == "PASS" ]]; then
        OVERALL_QC="$FILE_VALIDITY"
        return
    fi
    
    if [[ "$READ_QC" == "null" ]]; then
        OVERALL_QC="PREPROCESS MODULE FAILURE"
        return
    fi 
    
    if [[ "$READ_QC" == "FAIL" ]]; then
        OVERALL_QC="FAIL"
        return
    fi 

    OVERALL_QC=""

    if [[ "$ASSEMBLY_QC" == "null" ]]; then
        OVERALL_QC+="ASSEMBLY MODULE FAILURE;"
    fi
    
    if [[ "$MAPPING_QC" == "null" ]]; then
        OVERALL_QC+="MAPPING MODULE FAILURE;"
    fi

    if [[ "$TAXONOMY_QC" == "null" ]]; then
        OVERALL_QC+="TAXONOMY MODULE FAILURE;"
    fi

    if [[ ! "$OVERALL_QC" == "" ]]; then
        OVERALL_QC="${OVERALL_QC%;}"
        return
    fi

    if [[ "$READ_QC" == "PASS" ]] && [[ "$ASSEMBLY_QC" == "PASS" ]] && [[ "$MAPPING_QC" == "PASS" ]] && [[ "$TAXONOMY_QC" == "PASS" ]]; then
        OVERALL_QC="PASS"
    else
        OVERALL_QC="FAIL"
    fi
}

assign_overall_qc

echo \"Overall_QC\" > "$OVERALL_QC_REPORT"
echo \""$OVERALL_QC"\" >> "$OVERALL_QC_REPORT"


