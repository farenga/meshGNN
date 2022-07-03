Lanciare lo script main.m prima di cominciare e non spostare il file dalla
sua posizione orginale.

Nella cartella examples ci sono degli esempi per capire come funziona il
codice.

Script utili che andranno modificati durante il progetto (in ordine):
- mesh\build\grids2D.m
- mesh\agglomerate\aggl_GNN_fun.m
- numerical_experiments\quality\uniform_agglomeration.m
- numerical_experiments\quality\plot_aggl_grids.m
- numerical_experiments\quality\compute_quality.m
- numerical_experiments\quality\plot_quality.m
- numerical_experiments\multigrid\MG_test.m

I file di output come le griglie o le immagini vengono salvati nella
cartella output_files.

Per utilizzare il wrapper di python è necessario lanciare i comandi nel
file txt "comandi".

Per lanciare la funzione aggl_GNN_fun.m è necessario sostituire in 
'model directory' la directory al file model_res.pt o model_base.pt,
facendo attenzione alla riga 20 di inserire il file relativo.
