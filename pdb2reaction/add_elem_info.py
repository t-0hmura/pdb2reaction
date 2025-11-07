# pdb2reaction/add_elem_info.py

"""
PDB の元素列（columns 77–78）が無い/不完全なファイルに元素情報を付与します。
- Biopython で PDB を読み込み、atom.element を設定してから PDBIO で書き出し
- 原子名だけでなく残基名（タンパク/核酸/水/イオン）も併用して判定
- 出力パス未指定時は入力ファイルへ上書き保存
- 【新機能】入力に既に element 情報がある原子は、--overwrite を付けない限り「上書きせずにスキップ」

使い方:
  # 単体スクリプトとして
  python add_elem_info.py input.pdb [-o fixed.pdb] [--overwrite]

  # pdb2reaction サブコマンドとして
  pdb2reaction add_elem_info -i input.pdb [-o fixed.pdb] [--overwrite]
"""
from __future__ import annotations

import argparse
import collections
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Set

import click
from Bio.PDB import PDBParser, PDBIO

# -----------------------------
# 元素シンボル（IUPAC, 1–118）
# -----------------------------
ELEMENTS: Set[str] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db",
    "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Fl","Lv","Ts","Og"
}

# -----------------------------
# ご提示の辞書（そのまま使用）
# -----------------------------
AMINO_ACIDS: Dict[str, int] = {
    # --- Standard 20 (L) ---
    "ALA":  0, "ARG": +1, "ASN":  0, "ASP": -1, "CYS":  0,
    "GLU": -1, "GLN":  0, "GLY":  0, "HIS":  0, "ILE":  0,
    "LEU":  0, "LYS": +1, "MET":  0, "PHE":  0, "PRO":  0,
    "SER":  0, "THR":  0, "TRP":  0, "TYR":  0, "VAL":  0,

    # --- Canonical extras ---
    "SEC": -1,   # selenocysteine
    "PYL": +1,   # pyrrolysine

    # --- Protonation / tautomers (Amber/CHARMM style) ---
    "HIP": +1,   # fully protonated His
    "HID":  0,   # Nδ-protonated His
    "HIE":  0,   # Nε-protonated His
    "ASH":  0,   # neutral Asp
    "GLH":  0,   # neutral Glu
    "LYN":  0,   # neutral Lys

    # --- Phosphorylated residues ---
    "SEP": -2, "TPO": -2, "PTR": -2,

    # --- Cys family ---
    "CYX":  0,   # disulfide Cys
    "CSO":  0,   # Cys sulfenic acid
    "CSD": -1,   # Cys sulfinic acid
    "CSX":  0,   # generic Cys derivative
    "OCS": -1,   # cysteic acid

    # --- Lys variants / carboxylation ---
    "MLY": +1, "LLP": +1, "DLY": +1,
    "KCX": -1,   # Lysine Nz-Carboxylic Acid

    # --- D isomers (19 residues) ---
    "DAL":  0, "DAR": +1, "DSG": 0, "DAS": -1, "DCY": 0,
    "DGN":  0, "DGL": -1, "DHI": 0, "DIL":  0, "DLE": 0,
    "DLY": +1, "MED":  0, "DPN": 0, "DPR":  0, "DSN": 0,
    "DTH":  0, "DTR":  0, "DTY": 0, "DVA":  0,

    # --- Carboxylation / cyclization / others ---
    "CGU": -2,   # gamma-carboxy-glutamate
    "CGA": -1,   # carboxymethylated glutamate
    "PCA":  0,   # pyroglutamate
    "MSE":  0,   # selenomethionine
    "OMT":  0,   # methionine sulfone

    # --- Other modified residues possibly encountered ---
    "ASA": 0, "CIR": 0, "FOR": 0, "MVA": 0, "IIL": 0, "AIB": 0, "HTN": 0,
    "SAR": 0, "NMC": 0, "PFF": 0, "NFA": 0, "ALY": 0, "AZF": 0, "CNX": 0, "CYF": 0,

    # --- Hydroxyproline ---
    "HYP": 0,

    # --- All C-terminus ---
    "CALA": -1, "CARG":  0, "CASN": -1, "CASP": -2, "CCYS": -1,
    "CCYX": -1, "CGLN": -1, "CGLU": -2, "CGLY": -1, "CHID": -1,
    "CHIE": -1, "CHIP":  0, "CHYP": -1, "CILE": -1, "CLEU": -1,
    "CLYS":  0, "CMET": -1, "CPHE": -1, "CPRO": -1, "CSER": -1,
    "CTHR": -1, "CTRP": -1, "CTYR": -1, "CVAL": -1, "NHE": 0,
    "NME": 0,

    # --- All N-terminus ---
    "NALA": +1, "NARG": +2, "NASN": +1, "NASP":  0, "NCYS": +1,
    "NCYX": +1, "NGLN": +1, "NGLU":  0, "NGLY": +1, "NHID": +1,
    "NHIE": +1, "NHIP": +2, "NILE": +1, "NLEU": +1, "NLYS": +2,
    "NMET": +1, "NPHE": +1, "NPRO": +1, "NSER": +1, "NTHR": +1,
    "NTRP": +1, "NTYR": +1, "NVAL": +1, "ACE": 0,
}

ION: Dict[str, int] = {
    # +1
    "LI": +1, "NA": +1, "K": +1, "RB": +1, "CS": +1, "TL": +1, "AG": +1, "CU1": +1,
    "Ag": +1, "K+": +1, "Na+": +1, "NH4": +1, "H3O+": +1, "HE+": +1, "HZ+": +1, "Tl": +1,

    # +2
    "MG": +2, "CA": +2, "SR": +2, "BA": +2, "MN": +2, "FE2": +2, "CO": +2, "NI": +2,
    "CU": +2, "ZN": +2, "CD": +2, "HG": +2, "PB": +2, "Be": +2, "PD": +2, "PT": +2,
    "Sn": +2, "Ra": +2, "YB2": +2, "V2+": +2,

    # +3
    "FE": +3, "AU3": +3, "AL": +3, "GA": +3, "IN": +3,
    "CE": +3, "Ce": +3, "CR": +3, "Cr": +3, "Dy": +3, "EU": +3, "EU3": +3, "Er": +3,
    "GD3": +3, "LA": +3, "LU": +3, "Nd": +3, "PR": +3, "SM": +3, "Sm": +3, "TB": +3,
    "Tm": +3, "Y": +3, "Pu": +3,

    # +4
    "U4+": +4, "Th": +4, "Hf": +4, "Zr": +4,

    # -1
    "F": -1, "CL": -1, "BR": -1, "I": -1, "Cl-": -1, "IOD": -1,
}

# よく使う残基クラス
PROTEIN_RES = set(AMINO_ACIDS.keys())
NUCLEIC_RES = {
    # DNA/RNA（最低限）
    "DA","DT","DG","DC","DI",
    "A","U","G","C","I",
}
WATER_RES = {"HOH","WAT","H2O","DOD","TIP","TIP3","SOL"}

# -----------------------------
# ヘルパー: 文字列→元素シンボル正規化
# -----------------------------
_re_letters = re.compile(r"[A-Za-z]+")

def _normalize_symbol(s: str) -> Optional[str]:
    """英字以外を除去し、2文字優先で既知元素にマッチしたらシンボル（正しい大小）を返す。"""
    if not s:
        return None
    m = _re_letters.findall(s)
    if not m:
        return None
    letters = "".join(m)
    if len(letters) >= 2:
        cand2 = (letters[:2][0].upper() + letters[:2][1].lower())
        if cand2 in ELEMENTS:
            return cand2
    cand1 = letters[0].upper()
    if cand1 in ELEMENTS:
        return cand1
    # Deuterium を H として扱う（PDB では D を H と同等に扱うことが多い）
    if letters[0].upper() == "D":
        return "H"
    return None

def _symbol_from_resname(resname: str) -> Optional[str]:
    """イオン残基名（例: CA, FE2, Cl-, YB2, IOD）から元素シンボルを抽出。"""
    res = resname.strip()
    sym = _normalize_symbol(res)
    if sym is None and res.upper().startswith("IOD"):
        sym = "I"
    return sym

# -----------------------------
# 元素推定ロジック（残基名を使って曖昧さを解消）
# -----------------------------
def guess_element(atom_name: str, resname: str, is_het: bool) -> Optional[str]:
    """
    原子名 + 残基名 から元素を推定。
    優先順位:
      1) 残基がイオン: 残基名から元素（NH4/H3O+などの多原子種は原子名で H/N/O へ）
      2) タンパク/核酸/水: 慣例に従い H/C/N/O/S/P/Se を優先（CA=Carbon, HG=Hydrogen）
      3) その他リガンド: 原子名先頭で推定、ただし C* や P* パターンは Carbon/Phosphorus を優先
      4) それでも曖昧なら 2文字→1文字の順で元素候補にフォールバック
    """
    name_u = atom_name.strip().upper()
    res_u = resname.strip().upper()

    # 1) イオン残基（残基名→元素を強く優先）
    if res_u in {k.upper() for k in ION.keys()}:
        # 多原子イオン（NH4, H3O+, …）は原子名で個々判定（D* も H とみなす）
        if name_u.startswith(("H", "D")):
            return "H"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        # 単原子金属/ハロゲンなどは残基名から
        sym = _symbol_from_resname(res_u)
        if sym:
            return sym
        # 残基名が例外的でも、原子名が CL/BR/I/F ならそのハロゲン
        if name_u.startswith("CL"):
            return "Cl"
        if name_u.startswith("BR"):
            return "Br"
        if name_u.startswith("I"):
            return "I"
        if name_u.startswith("F"):
            return "F"

    # 2) ポリマー（タンパク/核酸）/水
    is_protein = res_u in PROTEIN_RES
    is_nucl = res_u in NUCLEIC_RES
    is_water = res_u in WATER_RES
    if is_protein or is_nucl or is_water:
        # 水は O と H のみ（D* も H とみなす）
        if is_water:
            if name_u.startswith(("H", "D")):
                return "H"
            return "O"

        # 水素（D* を含む）
        if name_u.startswith(("H", "D")):
            return "H"

        # セレノメチオニン/セレノシステイン等
        if name_u.startswith("SE"):
            return "Se"

        # P, N, O, S は素直に先頭文字で
        if name_u.startswith("P"):
            return "P"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        if name_u.startswith("S"):
            return "S"

        # Cα/側鎖（CA, CB, CG, CD, CE, CZ, CH* など）は Carbon
        if name_u.startswith("C"):
            return "C"

        # まれにハロゲンを持つ場合もあるため最後にフォールバック
        sym = _normalize_symbol(name_u)
        if sym:
            return sym

    # 3) 非ポリマー（リガンド/補欠分子族など）
    #    炭素/リン様ラベル（C*, P*）は Carbon/Phosphorus を優先（ただし CL は除外）
    if name_u.startswith("C") and not name_u.startswith("CL"):
        return "C"
    if name_u.startswith("P"):
        return "P"

    # 金属やハロゲンがそのまま原子名のことも多い（FE, ZN, MG, HG, CL, BR, I, F ...）
    sym = _normalize_symbol(name_u)
    if sym:
        return sym

    # 4) どうしても決め打てない場合は None（呼び出し側で警告）
    return None

# -----------------------------
# 入力 PDB の element 列が「元々」存在していたかを
# シリアル番号（columns 7–11）で判定するヘルパー
# -----------------------------
def scan_existing_elements_by_serial(pdb_path: str) -> Set[int]:
    """
    PDB 行を直接スキャンし、element フィールド（columns 77–78）が空ではない ATOM/HETATM
    のシリアル番号を返す。Biopython の自動推定に影響されない「原ファイルの実際の有無」を見る。
    """
    serials_with_elem: Set[int] = set()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 78:
                    # element 列が存在しない
                    continue
                serial_str = line[6:11].strip()
                elem_raw = line[76:78].strip()
                if not serial_str:
                    continue
                try:
                    serial = int(serial_str)
                except ValueError:
                    continue
                # 空でなければ「element 情報があった」とみなす（D などの同位体表記も保持したい）
                if elem_raw:
                    serials_with_elem.add(serial)
    except Exception:
        # 読めない場合は空集合（= すべて未設定扱い）
        pass
    return serials_with_elem

def _get_atom_serial(atom) -> Optional[int]:
    """Biopython Atom からシリアル番号を安全に取得。バージョン差異を吸収。"""
    sn = getattr(atom, "serial_number", None)
    if sn is None and hasattr(atom, "get_serial_number"):
        try:
            sn = atom.get_serial_number()
        except Exception:
            sn = None
    return sn

# -----------------------------
# メイン処理
# -----------------------------
def assign_elements(in_pdb: str, out_pdb: Optional[str], overwrite: bool = False) -> None:
    # 入力の element 列の「生有無」をスキャン
    existing_by_serial = scan_existing_elements_by_serial(in_pdb)

    parser = PDBParser(QUIET=True)
    structure_id = os.path.splitext(os.path.basename(in_pdb))[0]
    structure = parser.get_structure(structure_id, in_pdb)

    total = 0
    assigned_new = 0          # element フィールドが無かった原子に新規設定
    overwritten = 0           # element フィールドが元々あったが --overwrite により再設定
    kept_existing = 0         # element フィールドが元々あり、--overwrite 無しで保持
    unknown = []              # 推定不能（未変更）

    by_element = collections.Counter()

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()  # '' (空) なら標準、'W' は水、'H_' は HETATM
                is_het = (hetflag != "")
                resname = residue.get_resname()
                for atom in residue:
                    total += 1
                    name = atom.get_name()

                    serial = _get_atom_serial(atom)
                    had_element_in_input = (serial in existing_by_serial) if serial is not None else False

                    if had_element_in_input and not overwrite:
                        kept_existing += 1
                        continue  # 既存の element を尊重して一切変更しない

                    sym = guess_element(name, resname, is_het)
                    if sym is None:
                        unknown.append((model.id, chain.id, residue.id, resname, name, serial))
                        # 推定不可：既存値があるなら温存、なければ未設定のまま
                        continue

                    # Biopython は atom.element を見て列 77–78 を出力
                    prev = getattr(atom, "element", None)
                    atom.element = sym
                    by_element[sym] += 1
                    if had_element_in_input:
                        if prev != sym:
                            overwritten += 1
                    else:
                        assigned_new += 1

    io = PDBIO()
    io.set_structure(structure)
    out_path = out_pdb if out_pdb else in_pdb  # 指定なしは上書き
    io.save(out_path)

    # サマリ
    print(f"[OK] Wrote: {out_path}")
    print(f"  atoms total            : {total}")
    print(f"  newly assigned         : {assigned_new}")
    print(f"  kept existing (no ow)  : {kept_existing}")
    print(f"  overwritten (--overwrite): {overwritten}")
    if by_element:
        top = ", ".join(f"{k}:{v}" for k, v in by_element.most_common())
        print(f"  assigned breakdown     : {top}")
    if unknown:
        print(f"[WARN] Could not confidently assign {len(unknown)} atoms; left unchanged.")
        for (mid, chid, resid, resn, aname, serial) in unknown[:50]:
            if isinstance(resid, tuple):
                resseq = resid[1]
                icode = resid[2].strip()
            else:
                resseq, icode = "?", ""
            s_str = f" serial {serial}" if serial is not None else ""
            print(f"    model {mid} chain {chid} {resn} {resseq}{icode} : {aname}{s_str}")
        if len(unknown) > 50:
            print("    ... (truncated) ...")

def main():
    ap = argparse.ArgumentParser(
        description="Add/repair element columns (77–78) in a PDB using Biopython."
    )
    ap.add_argument("pdb", help="input PDB filepath")
    ap.add_argument("-o", "--out", help="output PDB filepath (omit to overwrite input)")
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="既に element 列がある原子も再推定して上書きする（指定が無い場合は既存を保持）",
    )
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        print(f"[ERR] Input not found: {args.pdb}", file=sys.stderr)
        sys.exit(1)

    try:
        assign_elements(args.pdb, args.out, overwrite=args.overwrite)
    except Exception as e:
        print(f"[ERR] Failed: {e}", file=sys.stderr)
        sys.exit(2)

# -----------------------------
# Click サブコマンド（pdb2reaction add_elem_info）
# -----------------------------
@click.command(
    help="Add/repair element columns (77–78) in a PDB using Biopython.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "in_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input PDB filepath",
)
@click.option(
    "-o", "--out",
    "out_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help="Output PDB filepath (omit to overwrite input)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="既に element 列がある原子も再推定して上書きする（指定が無い場合は既存を保持）",
)
def cli(in_pdb: Path, out_pdb: Optional[Path], overwrite: bool) -> None:
    """pdb2reaction のサブコマンド経由で実行するための Click ラッパー。"""
    try:
        assign_elements(str(in_pdb), (str(out_pdb) if out_pdb else None), overwrite=overwrite)
    except SystemExit as e:
        # argparse 互換の main() と行儀を合わせるため、SystemExit はそのまま伝播
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)

if __name__ == "__main__":
    main()
