"""STEP2-2 : GSE30529 (HG-U133A_2) — RMA 는 R 전용. R 산출 매트릭스 재사용/확인."""
import importlib.util, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
_note = pathlib.Path(__file__).with_name("_rma_note.py")
_spec = importlib.util.spec_from_file_location("_rma_note", _note)
rma = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rma)

if __name__ == "__main__":
    rma.load_labeled("GSE30529")   # 기대: 13041 genes x 22 (Control 12 / DKD 10)
