
rule mkdir:
    output:
        "{out}/mkdir_tmp.txt".format(out=config["dir"]["out"])
    run:
        out = config["dir"]["out"]
        directories = []
        for k, v in config["dir"].items():
            if k != "out":
                directories.append(out + '/' + v)
        directories = ' '.join(directories)
        shell("mkdir -p {directories}")
        shell("touch {output} ")
        shell("sleep 5")
