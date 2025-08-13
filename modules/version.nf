process get_version {
    publishDir "${params.output}/meta", mode: 'copy'

    output:
    path "version.txt"

    shell:
    '''
    set -euo pipefail

    output="version.txt"
    version="${params.version:-}"

    echo "version" > "${output}"

    if [[ -z "${version}" ]]; then
      # Try tag first
      v="$(git -C "${baseDir}" describe --tags 2>/dev/null || true)"
      # If no tag, use short commit hash
      if [[ -z "${v}" ]]; then
        v="$(git -C "${baseDir}" rev-parse --short HEAD 2>/dev/null || true)"
      fi
      # If not a git repo, fall back to "unknown"
      [[ -z "${v}" ]] && v="unknown"
      echo "${v}" >> "${output}"
    else
      echo "${version}" >> "${output}"
    fi
    '''
}