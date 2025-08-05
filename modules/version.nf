process get_version {
    input:

    output:
    path("${output}")

    script:
    output="version.txt"
    version=params.version
    """
    echo "version" > ${output}

    if [ -z ${version} ];
    then
        echo \$(git -C $baseDir describe --tags) >> ${output}
    else
        echo ${version} >> ${output}
    fi
    """
}
