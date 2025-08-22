process get_version {
    label 'farm_mid'
    
    output:
    path "version.txt"

    script:
    """
    set -euo pipefail

    echo "version" > version.txt

    if [ -z "${params.version}" ]; then
      v="\$(git -C "${workflow.projectDir}" describe --tags 2>/dev/null || true)"
      if [ -z "\$v" ]; then
        v="\$(git -C "${workflow.projectDir}" rev-parse --short HEAD 2>/dev/null || true)"
      fi
      [ -z "\$v" ] && v="unknown"
      echo "\$v" >> version.txt
    else
      echo "${params.version}" >> version.txt
    fi
    """
}
