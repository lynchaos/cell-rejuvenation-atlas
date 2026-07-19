"""Tests for the module-2 time-course loader (synthetic 10x h5 + name parsing)."""
import numpy as np
import pytest
import scipy.sparse as sp

h5py = pytest.importorskip("h5py")
pytest.importorskip("scanpy")

from src.module2_reprogramming_trajectory.preprocess import load_matrix, parse_sample_name


def test_parse_sample_name_variants():
    assert parse_sample_name("GSM2836276_D9-1-2i") == (9.0, "2i")
    assert parse_sample_name("GSM2836269_D2-2") == (2.0, "serum")
    assert parse_sample_name("GSM3195648_D0_Dox_C1_gene_bc_mat") == (0.0, "dox")
    assert parse_sample_name("GSM3195650_D0.5_Dox_C1_gene_bc_mat") == (0.5, "dox")
    day, _ = parse_sample_name("GSM3195800_DiPSC_C1_gene_bc_mat")
    assert np.isnan(day)


def _write_10x_h5(path, n_cells=6, n_genes=4, seed=0):
    """Minimal cellranger-v3-format h5 (matrix group with features)."""
    x = sp.csr_matrix(np.random.default_rng(seed).integers(0, 3, size=(n_cells, n_genes)))
    xt = x.T.tocsc()  # 10x stores genes x cells as CSC (indptr = n_cells + 1)
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.array([f"c{i}".encode() for i in range(n_cells)]))
        m.create_dataset("data", data=xt.data.astype(np.int32))
        m.create_dataset("indices", data=xt.indices.astype(np.int32))
        m.create_dataset("indptr", data=xt.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array([n_genes, n_cells]))
        feats = m.create_group("features")
        feats.create_dataset("id", data=np.array([f"g{i}".encode() for i in range(n_genes)]))
        feats.create_dataset("name", data=np.array([f"Gene{i}".encode() for i in range(n_genes)]))
        feats.create_dataset("feature_type", data=np.array([b"Gene Expression"] * n_genes))


def test_load_matrix_reads_dox_h5_time_course(tmp_path):
    _write_10x_h5(tmp_path / "GSM1_D0_Dox_C1_gene_bc_mat.h5", seed=1)
    _write_10x_h5(tmp_path / "GSM2_D0.5_Dox_C1_gene_bc_mat.h5", seed=2)
    _write_10x_h5(tmp_path / "GSM3_D18_Dox_C1_gene_bc_mat.h5", seed=3)
    _write_10x_h5(tmp_path / "GSM4_DiPSC_C1_gene_bc_mat.h5", seed=4)  # no day -> dropped
    adata = load_matrix(tmp_path)
    assert set(adata.obs["day"]) == {0.0, 0.5, 18.0}
    assert set(adata.obs["condition"]) == {"dox"}
    assert adata.n_obs == 18  # 3 kept samples x 6 cells
    assert list(adata.var_names) == ["Gene0", "Gene1", "Gene2", "Gene3"]


def test_load_matrix_falls_back_to_triplets(tmp_path):
    import gzip

    x = sp.csr_matrix(np.array([[1, 0], [0, 1], [2, 0]]))  # 3 cells x 2 genes
    xt = x.T.tocoo()  # genes x cells, as 10x mtx
    with gzip.open(tmp_path / "GSM9_D4-1.matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n%\n2 3 3\n")
        for i, j, v in zip(xt.row, xt.col, xt.data):
            f.write(f"{i + 1} {j + 1} {v}\n")
    with gzip.open(tmp_path / "GSM9_D4-1.barcodes.tsv.gz", "wt") as f:
        f.write("c0\nc1\nc2\n")
    with gzip.open(tmp_path / "GSM9_D4-1.genes.tsv.gz", "wt") as f:
        f.write("g0\tGene0\ng1\tGene1\n")
    adata = load_matrix(tmp_path)
    assert set(adata.obs["day"]) == {4.0}
    assert list(adata.var_names) == ["Gene0", "Gene1"]
