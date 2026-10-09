"""Peringatan bila data paper dan evaluator memakai faktor konversi ECOST yang berbeda.
paper_data/FAKTOR.txt ditulis oleh update_data.py setelah data dipindahkan ke faktor 1,7149."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "evaluator"))


def check():
    import ecost
    f = ecost.conversion_factor()
    data_new = os.path.exists(os.path.join(HERE, "paper_data", "FAKTOR.txt"))
    eval_new = abs(f - 1.7149) < 1e-3
    if data_new != eval_new:
        print(f"PERINGATAN: evaluator memakai faktor {f:.4f}, tetapi paper_data "
              f"{'sudah' if data_new else 'belum'} dipindahkan ke faktor 1,7149 (update_data.py). "
              "Nilai f2 hasil evaluasi ulang dan nilai f2 di berkas data tidak sebanding.", file=sys.stderr)
    return f
