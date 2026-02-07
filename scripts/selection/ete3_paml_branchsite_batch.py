from ete3 import EvolTree
import glob, os, csv
import shutil

# -----------------------------
ALIGN_DIR = "alignments/"
OUT_DIR = "results_model/"
os.makedirs(OUT_DIR, exist_ok=True)

TREE_FILE = "tree.nwk"  # 普通 Newick 树
FOREGROUND_LEAVES = ["Mitrastemon"]  # 可以标记多个叶子

# -----------------------------
for aln_file in glob.glob(os.path.join(ALIGN_DIR, "*.fasta")):
        gene = os.path.splitext(os.path.basename(aln_file))[0]
        workdir = os.path.join(OUT_DIR, gene)
        os.makedirs(workdir, exist_ok=True)

        try:
            # ---------- 初始化 EvolTree ----------
            tree = EvolTree(TREE_FILE)
            tree.link_to_alignment(aln_file)
            tree.workdir = workdir

            # ---------- 标记前景 ----------
            fg_ids = []
            for name in FOREGROUND_LEAVES:
                try:
                    fg_ids.append((tree & name).node_id)
                except Exception:
                    print(f"{gene}: Warning - leaf '{name}' not found in tree")
            if fg_ids:
                tree.mark_tree(fg_ids, ['#1'])
                print(f"{gene}: foreground branches marked")

            # ---------- 运行所有模型 ----------
            models = ["bsA1","bsA"]
            for model_name in models:
                try:
                    print(f"{gene}: running {model_name} ...")
                    # 不指定 out_file，Ete3 默认生成 'out'
                    tree.run_model(model_name)                    
                    print(f"{gene}: {model_name} finished")
                except Exception as e:
                    # 安全打印异常
                    if hasattr(e, 'args') and e.args:
                        err_msg = []
                        for a in e.args:
                            if isinstance(a, bytes):
                                err_msg.append(a.decode())
                            else:
                                err_msg.append(str(a))
                        err_msg = "; ".join(err_msg)
                    else:
                        err_msg = str(e)
                    print(f"{gene}: ERROR running {model_name}: {err_msg}")

            print(f"{gene}: DONE\n")

        except Exception as e:
            print(f"{gene}: ERROR unexpected: {e}")

