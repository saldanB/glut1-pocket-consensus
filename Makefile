.PHONY: all align consensus sensitivity iou verify clean
all: align consensus sensitivity verify
align:       ; python run/01_align.py
consensus:   ; python run/02_consensus.py
sensitivity: ; python run/03_sensitivity.py
quick:       ; python run/03_sensitivity.py --quick
iou:         ; python run/04_iou.py
verify:      ; python run/99_verify.py
clean:       ; rm -rf results
