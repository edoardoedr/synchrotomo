import yaml
import sys
from src.stp_core import exec_tiff2tdf
from src.stp_core import exec_phaseretrieval
from src.stp_core import exec_reconstruct


def run_tiff2tdf_from_config(config_file):
    """
    Esegue exec_tiff2tdf.py usando i parametri da un file di configurazione YAML.
    
    Parameters
    ----------
    config_file : str
        Percorso al file di configurazione YAML
        
    Example YAML structure:
    ------------------------
    from: 0
    to: -1
    input_path: "/path/to/input/images/"
    output_file: "/path/to/output.tdf"
    crop_top: 0
    crop_bottom: 0
    crop_left: 0
    crop_right: 0
    tomo_prefix: "tomo"
    flat_prefix: "flat"
    dark_prefix: "dark"
    projection_order: True
    privilege_sino: True
    compression: 5
    nr_threads: 4
    log_file: "/path/to/log.txt"
    """
    
    # Carica il file di configurazione YAML
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Controlla se tiff2tdf  abilitato
    if 'tiff2tdf' not in config or not config['tiff2tdf'].get('enabled', True):
        print("[1/3] TIFF to TDF disabilitato nella configurazione. Skip.")
        return config.get('tiff2tdf', {}).get('output_file', '')

    tiff2tdf_config = config['tiff2tdf']
    
    # Prepara gli argomenti come stringhe (come se fossero passati da riga di comando)
    argv = [
        str(tiff2tdf_config.get('from', 0)),
        str(tiff2tdf_config.get('to', -1)),
        tiff2tdf_config.get('input_path', ''),
        tiff2tdf_config.get('output_file', ''),
        str(tiff2tdf_config.get('crop_top', 0)),
        str(tiff2tdf_config.get('crop_bottom', 0)),
        str(tiff2tdf_config.get('crop_left', 0)),
        str(tiff2tdf_config.get('crop_right', 0)),
        tiff2tdf_config.get('tomo_prefix', 'tomo'),
        tiff2tdf_config.get('flat_prefix', 'flat'),
        tiff2tdf_config.get('dark_prefix', 'dark'),
        str(tiff2tdf_config.get('projection_order', True)),
        str(tiff2tdf_config.get('privilege_sino', True)),
        str(tiff2tdf_config.get('compression', 5)),
        str(tiff2tdf_config.get('nr_threads', 4)),
        tiff2tdf_config.get('log_file', 'conversion.log')
    ]
    
    # Esegue la conversione TIFF to TDF
    print("[1/3] Avvio conversione TIFF to TDF con configurazione da: {}".format(config_file))
    exec_tiff2tdf.main(argv)
    print("[1/3] Conversione TIFF to TDF completata!")
    
    return tiff2tdf_config.get('output_file', '')


def run_phaseretrieval_from_config(config_file, input_tdf_file):
    """
    Esegue exec_phaseretrieval.py usando i parametri da un file di configurazione YAML.
    
    Parameters
    ----------
    config_file : str
        Percorso al file di configurazione YAML
    input_tdf_file : str
        Percorso al file TDF di input (output dello step precedente)
        
    Phase Retrieval parameters in YAML:
    ------------------------------------
    phase_retrieval:
        enabled: True
        method: 0  # 0=TIE-Hom (Paganin 2002), 1=TIE-Hom (Paganin 2020), 2+=Projected CTF
        param1: 1.0  # regParam o beta
        param2: 2000.0  # threshold o delta
        energy: 22.0  # keV
        distance: 300.0  # mm
        pixsize: 2.2  # micron
        pad: True
        output_file: "/path/to/output_phrt.tdf"
        log_file: "/path/to/phrt_log.txt"
    """
    
    # Carica il file di configurazione YAML
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Controlla se il phase retrieval abilitato
    if 'phase_retrieval' not in config or not config['phase_retrieval'].get('enabled', False):
        print("[2/3] Phase retrieval disabilitato nella configurazione. Skip.")
        return input_tdf_file
    
    pr_config = config['phase_retrieval']
    
    # Prepara gli argomenti per phase retrieval
    argv = [
        str(config.get('from', 0)),
        str(config.get('to', -1)),
        input_tdf_file,  # input file (output dello step precedente)
        pr_config.get('output_file', input_tdf_file.replace('.tdf', '_phrt.tdf')),
        str(pr_config.get('method', 0)),
        str(pr_config.get('param1', 1.0)),
        str(pr_config.get('param2', 2000.0)),
        str(pr_config.get('energy', 22.0)),
        str(pr_config.get('distance', 300.0)),
        str(pr_config.get('pixsize', 2.2)),
        str(pr_config.get('pad', True)),
        str(config.get('nr_threads', 4)),
        pr_config.get('log_file', 'phaseretrieval.log')
    ]
    
    # Esegue il phase retrieval
    print("[2/3] Avvio phase retrieval...")
    exec_phaseretrieval.main(argv)
    print("[2/3] Phase retrieval completato!")
    
    return pr_config.get('output_file', input_tdf_file.replace('.tdf', '_phrt.tdf'))


def run_reconstruct_from_config(config_file, input_tdf_file):
    """
    Esegue exec_reconstruct.py usando i parametri da un file di configurazione YAML.
    
    Parameters
    ----------
    config_file : str
        Percorso al file di configurazione YAML
    input_tdf_file : str
        Percorso al file TDF di input (output dello step precedente)
    """
    
    # Carica il file di configurazione YAML
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Controlla se la ricostruzione  abilitata
    if 'reconstruct' not in config or not config['reconstruct'].get('enabled', False):
        print("[3/3] Ricostruzione disabilitata nella configurazione. Skip.")
        return
    
    rec_config = config['reconstruct']
    
    # Prepara gli argomenti per la ricostruzione (42 parametri!)
    argv = [
        str(rec_config.get('from', 0)),
        str(rec_config.get('to', -1)),
        input_tdf_file,  # input file TDF
        rec_config.get('output_path', './reconstructed/'),
        str(rec_config.get('angles', 3.1416)),  # angoli in radianti (180 gradi = pi)
        str(rec_config.get('offset', 0.0)),
        rec_config.get('param1', 'shepp-logan'),  # filtro per FBP
        str(rec_config.get('scale', 1.0)),
        str(rec_config.get('overpad', True)),
        str(rec_config.get('logtransform', True)),
        str(rec_config.get('circle', True)),
        rec_config.get('output_prefix', 'recon'),
        str(rec_config.get('preprocessing_required', False)),
        str(rec_config.get('flat_end', False)),
        str(rec_config.get('half_half', False)),
        str(rec_config.get('half_half_line', 0)),
        str(rec_config.get('ext_fov', False)),
        str(rec_config.get('ext_fov_rot_right', False)),
        str(rec_config.get('ext_fov_overlap', 0)),
        str(rec_config.get('ext_fov_normalize', False)),
        str(rec_config.get('ext_fov_average', False)),
        str(rec_config.get('ext_fov_proj', 0)),
        str(rec_config.get('norm_sx', 0)),
        str(rec_config.get('norm_dx', 0)),
        str(rec_config.get('skip_ringrem', True)),
        rec_config.get('ringrem', '-'),
        str(rec_config.get('zerone_mode', False)),
        str(rec_config.get('corr_offset', 0.0)),
        rec_config.get('method', 'FBP_CUDA'),
        str(rec_config.get('decim_factor', 1)),
        str(rec_config.get('downsc_factor', 1)),
        str(rec_config.get('postprocess_required', False)),
        rec_config.get('polarfilt_opt', '-'),
        rec_config.get('convert_opt', '-'),
        rec_config.get('crop_opt', '-'),
        str(rec_config.get('angles_projfrom', 0)),
        str(rec_config.get('angles_projto', 1799)),
        str(rec_config.get('rolling', False)),
        str(rec_config.get('roll_shift', 0)),
        str(rec_config.get('dynamic_ff', False)),
        str(config.get('tiff2tdf', {}).get('nr_threads', 4)),
        rec_config.get('log_file', 'reconstruction.log')
    ]
    
    # Esegue la ricostruzione
    print("[3/3] Avvio ricostruzione tomografica...")
    exec_reconstruct.main(argv)
    print("[3/3] Ricostruzione completata!")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
        
        # Step 1: Conversione TIFF to TDF
        tdf_file = run_tiff2tdf_from_config(config_file)
            
        # Step 2: Phase Retrieval (se abilitato)
        phrt_file = run_phaseretrieval_from_config(config_file, tdf_file)
        
        # Step 3: Ricostruzione (se abilitata)
        run_reconstruct_from_config(config_file, phrt_file if phrt_file else tdf_file)
        
        print("\n=== Pipeline completata con successo! ===")
    else:
        print("Uso: python main.py <config_file.yaml>")
        print("\nEsempio di file di configurazione YAML:")
        print("""
from: 0
to: -1
input_path: "/percorso/alle/immagini/"
output_file: "/percorso/output.tdf"
crop_top: 0
crop_bottom: 0
crop_left: 0
crop_right: 0
tomo_prefix: "tomo"
flat_prefix: "flat"
dark_prefix: "dark"
projection_order: True
privilege_sino: True
compression: 5
nr_threads: 4
log_file: "/percorso/log.txt"

# Phase Retrieval (opzionale)
phase_retrieval:
  enabled: True
  method: 0  # 0=TIE-Hom (Paganin 2002), 1=TIE-Hom (Paganin 2020), 2+=CTF
  param1: 1.0  # beta (per TIE-Hom) o threshold (per CTF)
  param2: 2000.0  # delta (per TIE-Hom) o regParam (per CTF)
  energy: 22.0  # keV
  distance: 300.0  # mm
  pixsize: 2.2  # micron
  pad: True
  output_file: "/percorso/output_phrt.tdf"
  log_file: "/percorso/phrt_log.txt"
        """)

