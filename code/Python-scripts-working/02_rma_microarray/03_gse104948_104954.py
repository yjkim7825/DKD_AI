"""STEP2-3 : GSE104948/104954 (brainarray ENTREZG, 2-플랫폼) — RMA 는 R 전용.
R 산출 매트릭스 재사용/확인."""
import importlib.util, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
_note = pathlib.Path(__file__).with_name("_rma_note.py")
_spec = importlib.util.spec_from_file_location("_rma_note", _note)
rma = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rma)

if __name__ == "__main__":
    rma.load_labeled("GSE104948")  # 기대: 12042 genes x 38 (Control 26 / DKD 12)
    rma.load_labeled("GSE104954")  # 기대: 12042 genes x 43 (Control 26 / DKD 17)
