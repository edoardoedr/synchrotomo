# synchrotomo

Repo con **solo** la pipeline di ricostruzione tomografica di dati CT da sincrotrone (TIFF → TDF → phase retrieval → slice ricostruite), derivata da STP-Core (Elettra/SYRMEP). Documentazione utente completa in `README.md`.

## Punti chiave
- Entry point: `python main_reconstruction.py config_reconstruction.yaml`. `main_reconstruction.py` traduce le sezioni YAML (`tiff2tdf`, `phase_retrieval`, `reconstruct`) negli `argv` posizionali attesi da `src_reconstruction/exec_*.py:main(argv)` — **l'ordine degli argomenti conta** (42 per `exec_reconstruct`).
- Codice **Python 2.7** (env conda `synchrotomo`, scipy 0.12, ASTRA 2.1.2, CUDA). Non introdurre sintassi solo-Python 3 (f-string, ecc.) nel codice della pipeline.
- Gli import interni usano il package `src_reconstruction.` → eseguire dalla radice della repo.
- Le opzioni stringa `ringrem`, `polarfilt_opt`, `convert_opt` vengono fatte `split(":")`: usare `"none:"` per disattivarle, non `"-"`. `crop_opt` è `"top:bottom:left:right"`.
- `packages_reconstruction/` contiene estensioni C/Cython (gridrec, MR-FBP, FISTA-TV) da compilare con `install_dependencies.sh`.
- `vendor/STP_1.6.3/` è il toolbox originale, solo riferimento: non modificarlo né importarlo.
- Dati (`*.tdf`, `*.tif`, log) sono esclusi da git via `.gitignore`.
