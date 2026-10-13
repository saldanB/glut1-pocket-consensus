.PHONY: all correlation figure1 align consensus sensitivity quick iou persistence membrane trajectories dbscan deeppocket phases table4 verify clean
all: correlation align consensus sensitivity iou persistence membrane trajectories dbscan deeppocket phases table4 figure1 verify
correlation:  ; python run/00_correlation.py
align:        ; python run/01_align.py
consensus:    ; python run/02_consensus.py
sensitivity:  ; python run/03_sensitivity.py
quick:        ; python run/03_sensitivity.py --quick
iou:          ; python run/04_iou.py
persistence:  ; python run/05_persistence.py
membrane:     ; python run/06_membrane.py
trajectories: ; python run/07_iou_trajectories.py 0.5 5.0 && python run/07_iou_trajectories.py 0.3 5.0 && python run/07_iou_trajectories.py 0.7 5.0
dbscan:       ; python run/08_dbscan_vs_hdbscan.py
deeppocket:   ; python run/09_deeppocket_configs.py
phases:       ; python run/10_phase_comparison.py
table4:       ; python run/11_table4.py
figure1:      ; python run/12_figure1.py
verify:       ; python run/99_verify.py
clean:        ; rm -rf results
