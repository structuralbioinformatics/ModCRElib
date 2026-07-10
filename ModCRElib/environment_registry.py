TOOLS = {
    "blastp": {
        "name": "BLAST+",
        "type": "prebuilt",
        "url": "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz",
        "binary_path": "ncbi-blast-2.16.0+/bin/blastp",
    },


    "clustalo": {
        "name": "Clustal Omega",
        "type": "binary",
        "url": "http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz",
        "binary_path": "clustal-omega-1.2.4/src/clustalo",
    },

    "hmmscan": {
        "name": "HMMER",
        "type": "source",
        "url": "http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz",
        "build_steps": [
            "./configure --prefix=$INSTALL_DIR",
            "make -j4",
        ],
        "binary_path": "hmmer-3.3.2/src/hmmscan",
        "install_type": "source",
        "entrypoint": "hmmscan"
    },

}



