"""
hprcPrepare.py

Downloads the per-chromosome HPRC Minigraph-Cactus graphs, converts each to
GFA, normalizes it into a DAG with this repository's preprocessing
(Preprocessing/mcToDag.py) and writes a report of what the normalization
changed, over all chromosomes.

Source: HPRC release 1.1, Minigraph-Cactus graph on GRCh38 (the same graph
as chromosome 22 so far), one file per chromosome:
  https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/
      minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr<N>.vg

=====================================================================
HOW TO RUN
=====================================================================

Run from the repository root (needs curl, vg and enough disk space; vg:
static binary from https://github.com/vgteam/vg/releases):

    # all chromosomes 1-22, X, Y (skips what is already done)
    nohup python3 -u scripts/hprcPrepare.py > hprc_prepare.log 2>&1 &

    # some chromosomes, vg not on the PATH, 8 threads for vg convert
    python3 -u scripts/hprcPrepare.py -c 21 22 X --vg ~/tools/vg --threads 8

    # only list what would be done (sizes, free disk space)
    python3 scripts/hprcPrepare.py --dry-run

    # only (re)write the report from the existing stats.json files
    python3 scripts/hprcPrepare.py --report-only

Options:
    -c / --chromosomes  1 ... 22 X Y (also chr1 / c1)      (default: all 24)
    --variant           default (chr<N>.vg) or d9 (chr<N>.d9.vg)  (default: default)
    --ref               reference sample in the walks      (default: GRCh38)
    --vg                vg program                         (default: vg)
    --threads           threads for vg convert             (default: 4)
    --connections       parallel connections per download  (default: 64;
                        the HPRC bucket allows about 15 KB/s per connection,
                        so 64 connections give about 1 MB/s; 1 = one curl)
    --python            interpreter for mcToDag.py         (default: pypy3 if
                        installed, else python3; pypy3 is several times faster)
    --keep-vg           keep the downloaded .vg file (default: deleted after
                        the conversion; the GFA is kept)
    --force             redo every step, even if its result exists
    --dry-run           only print the plan
    --report-only       only write the report

=====================================================================
WHAT IT DOES (per chromosome, in input/HPRC/chromosome_<N>/Graphs/)
=====================================================================

  1. download   chr<N>.vg with curl in --connections byte ranges at once
                (into chr<N>.vg.parts/, joined at the end; resumable: an
                interrupted download continues), size checked against the
                server
  2. convert    vg convert -f -> chromosome_<N>.gfa.gz (written to a .part
                file first, renamed when complete)
  3. DAG        Preprocessing/mcToDag.py <gfa> --ref GRCh38 --pipeline-dir
                <Graphs> --name chromosome_<N>, the same call as
                scripts/pipelineRunner.py -d HPRC: numeric_edgelist.txt,
                reduced_graph.gfa, unique_elements.txt, ordered_haplotypes.txt,
                haplotype_removal_report.txt, chromosome_<N>.stats.json,
                chromosome_<N>.nodemap.tsv (with reference coordinates),
                chromosome_<N>.dag.gfa, chromosome_<N>.batches.tsv
  A step whose result exists is skipped (unless --force), so the script can
  be started again after an interruption. A chromosome that fails is
  reported and the script continues with the next one.

  Before a download the free disk space is checked: at least 4 times the
  file size (the GFA is larger than the .vg, and mcToDag writes more files).

=====================================================================
REPORT
=====================================================================

  output/HPRC/hprc_dag_report.csv / .tex, one row per chromosome (from
  chromosome_<N>.stats.json): genomes, nodes / links / bp of the graph,
  links removed (inversion links, cycle-breaking back edges, links dropped
  with small components; all and in percent), nodes removed (and bp),
  walk transitions broken, the DAG's nodes / edges / sources / sinks,
  acyclic. This says how faithful each DAG is to its graph.
"""
import argparse
import csv
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
BASE_URL = ("https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/"
            "minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/")
ALL_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def chrom_name(c):
    """'1', 'chr1', 'c1', 'chromosome_1' -> '1'"""
    c = str(c)
    for prefix in ("chromosome_", "chr", "c"):
        if c.lower().startswith(prefix):
            c = c[len(prefix):]
            break
    return c.upper() if c.lower() in ("x", "y") else c


def remote_size(url):
    out = subprocess.run(["curl", "-sIL", "--max-time", "60", url], capture_output=True, text=True)
    size = None
    for line in out.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            size = int(line.split(":", 1)[1].strip())
    return size


def human(n):
    return f"{n / 1e9:.2f} GB" if n is not None else "?"


def download(url, target, size, connections=1):
    """Downloads url to target; with connections > 1 and a known size in that
    many byte ranges at once (the HPRC bucket limits every connection to
    about 15 KB/s, the speed grows with the number of connections)."""
    if connections > 1 and size is not None:
        return download_parallel(url, target, size, connections)
    part = target.with_name(target.name + ".part")
    for attempt in range(1, 6):
        r = subprocess.run(["curl", "-sS", "-L", "--fail", "--retry", "5", "-C", "-", "-o", str(part), url])
        have = part.stat().st_size if part.exists() else 0
        if r.returncode == 0 and (size is None or have == size):
            part.rename(target)
            return
        log(f"  download incomplete ({human(have)} of {human(size)}, curl exit {r.returncode}), retry {attempt}")
        time.sleep(10)
    raise RuntimeError(f"download of {url} failed")


def download_parallel(url, target, size, connections):
    """Range i of the file goes to <target>.parts/<i>. Each range is resumable:
    a curl only fetches what its file still lacks and appends it (a stalled
    connection, below 1 KB/s for 2 minutes, is ended and started again in the
    next round). When all ranges are complete they are joined in order into
    target and the size is checked."""
    parts_dir = target.with_name(target.name + ".parts")
    parts_dir.mkdir(exist_ok=True)
    single = target.with_name(target.name + ".part")
    if single.exists():                       # left over from a single-connection download
        single.unlink()
    chunk = -(-size // connections)
    ranges = [(i * chunk, min(size, (i + 1) * chunk) - 1) for i in range(connections) if i * chunk < size]
    files = [parts_dir / f"{i:03d}" for i in range(len(ranges))]

    def have(i):
        return files[i].stat().st_size if files[i].exists() else 0

    def total():
        return sum(have(i) for i in range(len(ranges)))

    for round_ in range(1, 31):
        for i, (lo, hi) in enumerate(ranges):
            if have(i) > hi - lo + 1:             # should not happen: start the range again
                files[i].unlink()
        todo = [i for i, (lo, hi) in enumerate(ranges) if have(i) < hi - lo + 1]
        if not todo:
            break
        log(f"  round {round_}: {len(todo)} of {len(ranges)} ranges to fetch, "
            f"{human(total())} of {human(size)} done")
        procs = []
        for i in todo:
            lo, hi = ranges[i]
            out = open(files[i], "ab")
            p = subprocess.Popen(["curl", "-sS", "-L", "--fail", "--speed-limit", "1024", "--speed-time", "120",
                                  "-r", f"{lo + have(i)}-{hi}", url], stdout=out)
            procs.append((p, out))
        start, start_bytes, last = time.time(), total(), time.time()
        while any(p.poll() is None for p, _ in procs):
            time.sleep(5)
            if time.time() - last >= 300:          # progress every 5 minutes
                last = time.time()
                done = total()
                rate = (done - start_bytes) / max(time.time() - start, 1)
                eta = (size - done) / rate / 60 if rate > 0 else float("inf")
                log(f"  {human(done)} of {human(size)} ({100 * done / size:.1f}%), "
                    f"{rate / 1024:.0f} KB/s, about {eta:.0f} min left")
        for p, out in procs:
            out.close()
        failed = sum(1 for p, _ in procs if p.returncode != 0)
        if failed:
            log(f"  {failed} connections ended early (curl errors or stalled), they continue in the next round")
            time.sleep(10)
    else:
        raise RuntimeError(f"download of {url} incomplete after 30 rounds ({human(total())} of {human(size)})")

    # join in order; each range is deleted once appended (little extra disk space)
    joined = target.with_name(target.name + ".joining")
    with open(joined, "wb") as out:
        for f in files:
            with open(f, "rb") as inp:
                shutil.copyfileobj(inp, out, 16 * 1024 * 1024)
            f.unlink()
    got = joined.stat().st_size
    if got != size:
        joined.unlink()
        raise RuntimeError(f"joined download has {got} bytes, expected {size}")
    joined.rename(target)
    parts_dir.rmdir()


def convert(vg, vg_file, gfa_gz, threads):
    part = gfa_gz.with_name(gfa_gz.name + ".part")
    with open(part, "wb") as out:
        p1 = subprocess.Popen([vg, "convert", "-f", "-t", str(threads), str(vg_file)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip", "-c"], stdin=p1.stdout, stdout=out)
        p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
    if rc1 != 0 or rc2 != 0:
        raise RuntimeError(f"vg convert (exit {rc1}) / gzip (exit {rc2}) failed for {vg_file}")
    part.rename(gfa_gz)


def prepare(c, args, python):
    name = f"chromosome_{c}"
    graphs = REPO / "input" / "HPRC" / name / "Graphs"
    vg_name = f"chr{c}.vg" if args.variant == "default" else f"chr{c}.d9.vg"
    url = BASE_URL + vg_name
    vg_file = graphs / vg_name
    gfa_gz = graphs / f"{name}.gfa.gz"
    gfa_plain = graphs / f"{name}.gfa"
    stats = graphs / f"{name}.stats.json"

    have_gfa = gfa_gz.exists() or gfa_plain.exists()
    gfa = gfa_gz if gfa_gz.exists() or not gfa_plain.exists() else gfa_plain
    have_dag = stats.exists() and (graphs / "numeric_edgelist.txt").exists() and \
        (graphs / f"{name}.nodemap.tsv").exists() and have_gfa and \
        stats.stat().st_mtime >= gfa.stat().st_mtime

    if args.dry_run:
        size = None if have_gfa and not args.force else remote_size(url)
        steps = []
        if args.force or not have_gfa:
            steps += [f"download {vg_name} ({human(size)})", "convert to GFA"]
        if args.force or not have_dag:
            steps.append("DAG (mcToDag.py)")
        log(f"chr{c}: " + (", ".join(steps) if steps else "done, nothing to do"))
        return "planned"

    graphs.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1 + 2: download and convert
    if args.force or not have_gfa:
        if not vg_file.exists() or args.force:
            size = remote_size(url)
            if size is None:
                raise RuntimeError(f"{url} not found")
            free = shutil.disk_usage(graphs).free
            if free < 4 * size:
                raise RuntimeError(f"not enough disk space for chr{c}: {human(free)} free, "
                                   f"need about {human(4 * size)}")
            log(f"chr{c}: downloading {vg_name} ({human(size)})")
            download(url, vg_file, size, args.connections)
        log(f"chr{c}: converting to GFA")
        convert(args.vg, vg_file, gfa_gz, args.threads)
        if not args.keep_vg:
            vg_file.unlink()
        gfa = gfa_gz
        log(f"chr{c}: {gfa.name} {human(gfa.stat().st_size)}")
    else:
        log(f"chr{c}: GFA exists ({gfa.name}), skipping download")

    # 3: DAG
    if args.force or not have_dag:
        log(f"chr{c}: normalizing to a DAG ({python} mcToDag.py)")
        cmd = [python, str(REPO / "Preprocessing" / "mcToDag.py"), str(gfa), "--ref", args.ref,
               "--pipeline-dir", str(graphs), "--name", name]
        with open(graphs / "mcToDag.log", "w") as lf:
            r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            raise RuntimeError(f"mcToDag.py failed (exit {r.returncode}), see {graphs / 'mcToDag.log'}")
    else:
        log(f"chr{c}: DAG exists, skipping")
    log(f"chr{c}: done in {time.time() - t0:.0f} s")
    return "ok"


REPORT_COLUMNS = [
    ("chromosome", None), ("genomes", ("walks", "batches")),
    ("nodes", ("input", "nodes")), ("links", ("input", "links")), ("bp", ("input", "bp")),
    ("inversion_links", ("links_removed", "inversion_links")), ("pct_inversion", ("links_removed", "pct_inversion")),
    ("back_edges", ("links_removed", "cycle_breaking_back_edges")),
    ("dropped_small_components", ("links_removed", "dropped_with_small_components")),
    ("links_removed", ("links_removed", "total_removed")), ("pct_links_removed", ("links_removed", "pct_total_removed")),
    ("nodes_removed", ("nodes_removed", "nodes")), ("pct_nodes_removed", ("nodes_removed", "pct_nodes")),
    ("bp_removed", ("nodes_removed", "bp")), ("pct_bp_removed", ("nodes_removed", "pct_bp")),
    ("broken_transitions", ("walks", "broken_transitions")), ("pct_broken", ("walks", "pct_broken")),
    ("dag_nodes", ("output", "nodes")), ("dag_edges", ("output", "edges")),
    ("sources", ("output", "sources")), ("sinks", ("output", "sinks")), ("acyclic", ("output", "acyclic")),
]


def write_report(chromosomes):
    rows = []
    for c in chromosomes:
        f = REPO / "input" / "HPRC" / f"chromosome_{c}" / "Graphs" / f"chromosome_{c}.stats.json"
        if not f.exists():
            continue
        st = json.load(open(f))
        row = {}
        for col, key in REPORT_COLUMNS:
            row[col] = f"chr{c}" if key is None else st.get(key[0], {}).get(key[1], "-")
        rows.append(row)
    if not rows:
        log("report: no stats.json found")
        return
    out = REPO / "output" / "HPRC"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "hprc_dag_report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[c for c, _ in REPORT_COLUMNS])
        w.writeheader()
        w.writerows(rows)
    # LaTeX: the columns a reader needs to judge the DAGs
    tex_cols = [("chromosome", "Chr."), ("genomes", "Genomes"), ("nodes", "Nodes"), ("links", "Links"),
                ("pct_inversion", "Inversion links (\\%)"), ("back_edges", "Back edges"),
                ("pct_links_removed", "Links removed (\\%)"), ("pct_nodes_removed", "Nodes removed (\\%)"),
                ("pct_bp_removed", "bp removed (\\%)"), ("pct_broken", "Walk transitions broken (\\%)"),
                ("dag_nodes", "DAG nodes"), ("dag_edges", "DAG edges")]
    lines = ["\\begin{table}[ht]", "\\centering", "\\begin{tabular}{|" + "c|" * len(tex_cols) + "}", "\\hline",
             " & ".join(f"\\textbf{{{h}}}" for _, h in tex_cols) + " \\\\", "\\hline"]
    for r in rows:
        lines.append(" & ".join(str(r[k]) for k, _ in tex_cols) + " \\\\")
        lines.append("\\hline")
    lines += ["\\end{tabular}",
              "\\caption{HPRC Minigraph-Cactus graphs (release 1.1, GRCh38) normalized to DAGs: what the "
              "normalization removed (inversion links, cycle-breaking back edges, small components) and how "
              "many haplotype walk transitions it breaks.}",
              "\\label{tab:hprc_dag_report}", "\\end{table}"]
    with open(out / "hprc_dag_report.tex", "w") as f:
        f.write("\n".join(lines) + "\n")
    log(f"report: {len(rows)} chromosomes -> {out / 'hprc_dag_report.csv'} and .tex")
    for r in rows:
        log(f"  {r['chromosome']:6s} genomes {r['genomes']}, nodes {r['nodes']}, links removed "
            f"{r['pct_links_removed']}%, nodes removed {r['pct_nodes_removed']}%, "
            f"transitions broken {r['pct_broken']}%, acyclic {r['acyclic']}")


def main():
    ap = argparse.ArgumentParser(description="Download HPRC Minigraph-Cactus chromosome graphs and make DAGs")
    ap.add_argument("-c", "--chromosomes", nargs="+", default=ALL_CHROMOSOMES)
    ap.add_argument("--variant", choices=["default", "d9"], default="default")
    ap.add_argument("--ref", default="GRCh38")
    ap.add_argument("--vg", default="vg")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--connections", type=int, default=64)
    ap.add_argument("--python", default=None)
    ap.add_argument("--keep-vg", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    args = ap.parse_args()
    chromosomes = [chrom_name(c) for c in args.chromosomes]
    python = args.python or (shutil.which("pypy3") or "python3")

    if not args.report_only:
        missing = [t for t in ("curl", "gzip") if shutil.which(t) is None]
        if not args.dry_run and shutil.which(args.vg) is None and not os.path.exists(args.vg):
            missing.append(f"vg ({args.vg})")
        if missing:
            sys.exit("missing programs: " + ", ".join(missing))
        log(f"chromosomes: {' '.join(chromosomes)}; mcToDag with {python}; "
            f"free disk space: {human(shutil.disk_usage(REPO).free)}")
        status = {}
        for c in chromosomes:
            try:
                status[c] = prepare(c, args, python)
            except Exception as e:
                status[c] = f"FAILED: {e}"
                log(f"chr{c}: FAILED: {e}")
        print("\n===================\nSUMMARY\n===================")
        for c in chromosomes:
            print(f"chr{c:5s} {status[c]}")
        if args.dry_run:
            return
    write_report(chromosomes)
    if not args.report_only and any(s.startswith("FAILED") for s in status.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
