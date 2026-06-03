# ===================================================
# Makefile для GROMACS симуляции тетрамера Abeta39_G25D
# ===================================================

.PHONY: all setup run_md process protein_only clean

# Основные цели:
#   make setup         — Полная подготовка системы (до создания md.tpr)
#   make run_md        — Запуск production run в screen
#   make process       — Основная обработка траектории (PBC, nojump, centering)
#   make protein_only  — Удаление растворителя → protein_only.xtc/tpr
#   make all           — Выполнить setup + run_md

all: setup run_md

# ====================== ПОДГОТОВКА СИСТЕМЫ ======================
setup: md.tpr

md.tpr: npt.gro
	gmx grompp -f ../../mdps/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2

npt.gro: nvt.gro
	gmx grompp -f ../../mdps/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
	gmx mdrun -deffnm npt

nvt.gro: em.gro
	gmx grompp -f ../../mdps/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
	gmx mdrun -deffnm nvt

em.gro: solvated_ions.pdb
	gmx grompp -f ../../mdps/minim.mdp -c solvated_ions.pdb -p topol.top -o em.tpr -maxwarn 2
	gmx mdrun -v -deffnm em

solvated_ions.pdb: ions.tpr
	echo "SOL" | gmx genion -s ions.tpr -o solvated_ions.pdb -p topol.top -pname NA -nname CL -neutral -conc 0.15

ions.tpr: solvated.pdb
	gmx grompp -f ../../mdps/ions.mdp -c solvated.pdb -p topol.top -o ions.tpr -maxwarn 2

solvated.pdb: boxed.pdb
	gmx solvate -cp boxed.pdb -cs spc216.gro -o solvated.pdb -p topol.top

boxed.pdb: processed.pdb
	gmx editconf -f processed.pdb -o boxed.pdb -c -d 1.2 -bt dodecahedron

processed.pdb:
	ln -sfn ../../charmm36-feb2026_cgenff-5.0.ff .
	gmx pdb2gmx -f ../tetramer_G25D.pdb -o processed.pdb -p topol.top -water tip3p -ignh -ff charmm36-feb2026_cgenff-5.0

# ====================== PRODUCTION MD ======================
run_md:
	screen -S md gmx mdrun -deffnm md -v

# ====================== ПРОЦЕССИНГ ТРАЕКТОРИИ ======================
# Центрирование по белку, nojump
process: md_clean.xtc

md_clean.xtc: step2.xtc md.tpr
	echo -e "1\n0" | gmx trjconv -f step2.xtc -s md.tpr -pbc mol -center -ur compact -o md_clean.xtc

step2.xtc: step1.xtc
	echo 0 | gmx trjconv -f step1.xtc -s md.tpr -pbc nojump -o step2.xtc

step1.xtc: md.xtc
	echo 0 | gmx trjconv -f md.xtc -s md.tpr -pbc whole -o step1.xtc

# ====================== ДОПОЛНИТЕЛЬНАЯ ОБРАБОТКА ======================
# Убрать растворитель
protein_only: protein_only.xtc protein_only.tpr

protein_only.xtc: md_clean.xtc
	echo -e "1\n1" | gmx trjconv -f md_clean.xtc -s md.tpr -o protein_only.xtc -pbc mol -center -ur compact

protein_only.tpr: md.tpr
	echo "1" | gmx convert-tpr -s md.tpr -o protein_only.tpr

# ====================== ОЧИСТКА ======================
clean:
	rm -f \#* *# *.tpr *.gro *.xtc *.cpt *.log *.edr *.xvg *.pdb *.itp *.trr *.top  mdout.mdp
