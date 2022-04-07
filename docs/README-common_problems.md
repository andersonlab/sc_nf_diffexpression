# Common problems

* If running on a virtual machine, you may need to set `export QT_QPA_PLATFORM="offscreen"` for `scanpy` to work properly. Issue described [here](https://github.com/ipython/ipython/issues/10627).
* We found Nextflow was writing some output into the `${HOME}` directory. This resulted in a Java error as soon as a Nextflow command was executed. Based on file sizes within `${HOME}`, it seemed like the output was being written within the conda environment (following `du -h | sort -V -k 1`). By deleting and re-installing the conda environment, the problem was solved. The below flags may help prevent this from the future. In addition, setting the flag `export JAVA_OPTIONS=-Djava.io.tmpdir=/path/with/enough/space/` may also help.

```bash
# To be run from the execution dir, before the above nextflow command
# If you are running this on a cluster, make sure you log into an interactive
# session with >25Gb of RAM.
export NXF_OPTS="-Xms25G -Xmx25G"
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/.nexflow_work"
export NXF_TEMP="${NXF_HOME}/.nexflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/.nexflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/.nexflow_singularity"
```
