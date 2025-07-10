#!/usr/bin/env bash

# Required inputs:
# - $BRACKEN_REPORT
# - $QC_SAGALACTIAE_PERCENTAGE
# - $QC_TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE
# - $TAXONOMY_QC_REPORT

# Get combined S. agalactiae abundance (% = fraction * 100)
PERCENTAGE=$(awk -F"\t" '
    $3 == "S" && ($1 == "Streptococcus agalactiae" || $1 == "Streptococcus sp. '\''group B'\''") { sum += $7 }
    END { printf "%.2f", sum * 100 }
' "$BRACKEN_REPORT")

# Get top non-agalactiae species (by highest fraction of total reads)
TOP_NON_AGALACTIAE_RECORD=$(awk -F"\t" '
    $3 == "S" && $1 != "Streptococcus agalactiae" && $1 != "Streptococcus sp. '\''group B'\''" { print }
' "$BRACKEN_REPORT" | sort -t $'\t' -k7,7nr | head -n 1)


if [[ -n "$TOP_NON_AGALACTIAE_RECORD" ]]; then
    TOP_NON_AGALACTIAE_SPECIES=$(echo "$TOP_NON_AGALACTIAE_RECORD" | cut -f1)
    RAW_ABUNDANCE=$(echo "$TOP_NON_AGALACTIAE_RECORD" | cut -f7)
    TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE=$(printf "%.3f" "$(echo "$RAW_ABUNDANCE * 100" | bc -l)")
else
    TOP_NON_AGALACTIAE_SPECIES="NONE_DETECTED"
    TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE="0.00"
fi

# Ensure main abundance is set if missing
PERCENTAGE=${PERCENTAGE:-0.00}

# QC logic
if [[ "$(echo "$PERCENTAGE >= $QC_SAGALACTIAE_PERCENTAGE" | bc -l)" == 1 ]] && \
   [[ "$(echo "$TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE <= $QC_TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE" | bc -l)" == 1 ]]; then
    TAXONOMY_QC="PASS"
else
    TAXONOMY_QC="FAIL"
fi

# Write CSV output
echo "\"Taxonomy_QC\",\"S.agalactiae_%\",\"Top_Non-agalactiae_Species\",\"Top_Non-agalactiae_Species_%\"" > "$TAXONOMY_QC_REPORT"
echo "\"$TAXONOMY_QC\",\"$PERCENTAGE\",\"$TOP_NON_AGALACTIAE_SPECIES\",\"$TOP_NON_AGALACTIAE_SPECIES_PERCENTAGE\"" >> "$TAXONOMY_QC_REPORT"
