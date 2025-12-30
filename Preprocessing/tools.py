import pathlib

def get_tools(repo_root):
    """
    Central registry of external tool paths.
    """
    return {
        "vg": str(repo_root / "Pangenome_Construction/vgTool/vg/bin/./vg"),
        "minichain": str(repo_root / "Pangenome_Construction/minichain/./minichain"),
        "minigraph": str(repo_root / "Pangenome_Construction/minigraph/./minigraph"),
        "mash": str(repo_root / "Pangenome_Construction/Mash/./mash"),
    }
