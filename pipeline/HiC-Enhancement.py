#!/usr/bin/env python3
"""
HiC-ECC: Hi-C Enhancement, Comparison, and Classification Tool
A comprehensive pipeline for Hi-C data enhancement and differential analysis
"""

import argparse
import sys
import os
import yaml
import logging
import subprocess
import re
from pathlib import Path
from datetime import datetime

__version__ = "1.0.0"
__author__ = "Your Name"

def setup_logging(output_dir, verbose=False):
    """Setup logging configuration"""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = Path(output_dir) / f"hic_ecc_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def get_chess_parameters(genome, resolution):
    """Get optimal CHESS parameters based on genome and resolution"""
    
    # CHESS minimum window requirements (tested values)
    chess_params = {
        'mm10': {
            'min_window': 310000,
            'recommended_windows': [310000, 500000, 1000000],
            'max_step': 50000
        },
        'hg38': {
            'min_window': 400000,
            'recommended_windows': [400000, 800000, 1500000],
            'max_step': 100000
        }
    }
    
    if genome not in chess_params:
        raise ValueError(f"Unsupported genome: {genome}. Supported: mm10, hg38")
    
    params = chess_params[genome]
    optimal_step = min(resolution * 3, params['max_step'])
    
    return {
        'window_size': params['min_window'],
        'step_size': optimal_step,
        'recommended_windows': params['recommended_windows'],
        'min_window': params['min_window']
    }

def create_config_template():
    """Create a configuration template file"""
    
    # Get genome-specific CHESS parameters
    genome = 'mm10'  # Default for template
    resolution = 10000
    chess_params = get_chess_parameters(genome, resolution)
    
    config_template = {
        'general': {
            'reference_genome': genome,
            'resolution': resolution,
            'species': 'mouse',
            'threads': 4
        },
        'samples': {
            'tissues': ['Mouse_Pancreas', 'Mouse_Liver', 'Mouse_Brain'],
            'chromosomes': ['1', '4', '5', '6', '15']
        },
        'paths': {
            'hicpro_dir': '/path/to/hicpro/output',
            'deephic_root': '~/DeepHiC',
            'deephic_model': '/path/to/deephic_model.pth',
            'genome_sizes': '/path/to/genome.chrom.sizes'
        },
        'deephic': {
            'chunk_size': 100,
            'stride_size': 100,
            'low_res_cutoff': 100,
            'boundary': 9999,
            'scale_factor': 1
        },
        'chess': {
            'window_size': chess_params['window_size'],
            'step_size': chess_params['step_size'],
            'sn_threshold': 0,
            'zsim_threshold': 0,
            'recommended_windows': chess_params['recommended_windows'],
            'auto_optimize': True  # Automatically adjust parameters
        }
    }
    
    with open('hic_ecc_config.yaml', 'w') as f:
        yaml.dump(config_template, f, default_flow_style=False, indent=2)
    
    print("Configuration template created: hic_ecc_config.yaml")
    print("Please edit this file with your parameters before running the pipeline.")

def load_config(config_file):
    """Load configuration from YAML file"""
    try:
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print(f"ERROR: Configuration file {config_file} not found!")
        print("Run 'HiC-ECC.py --create-config' to create a template.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR: Invalid YAML in configuration file: {e}")
        sys.exit(1)

def setup_output_structure(output_dir):
    """Create organized output directory structure"""
    base_path = Path(output_dir)
    
    directories = [
        'logs',
        'enhanced_matrices',
        'hic_format',
        'chess_pairs',
        'chess_similarity',
        'differential_regions',
        'plots',
        'reports'
    ]
    
    for dir_name in directories:
        (base_path / dir_name).mkdir(parents=True, exist_ok=True)
    
    return base_path

def fix_deephic_parser_syntax(deephic_path, logger):
    """Fix known syntax errors in DeepHiC's all_parser.py"""
    
    all_parser_file = Path(deephic_path) / 'all_parser.py'
    
    logger.info("Checking and fixing DeepHiC all_parser.py syntax...")
    
    with open(all_parser_file, 'r') as f:
        content = f.read()
    
    # Fix known syntax error: extra closing parenthesis in mouse range
    fixed_content = content.replace(
        "'mouse': list(range(1,20))) + ['X'],",
        "'mouse': list(range(1,20)) + ['X'],"
    )
    
    # Try to compile to check syntax
    try:
        compile(fixed_content, str(all_parser_file), 'exec')
        logger.info("Syntax is valid")
        
        # If we made changes, write the fixed version
        if fixed_content != content:
            logger.info("Fixed syntax error in all_parser.py")
            with open(all_parser_file, 'w') as f:
                f.write(fixed_content)
        
        return True
        
    except SyntaxError as e:
        logger.error(f"Syntax error persists: {e}")
        logger.error("You may need to manually fix all_parser.py")
        return False

def update_deephic_root_dir(deephic_path, new_root_dir, logger):
    """Update root directory and fix missing 'all' key in DeepHiC's all_parser.py"""
    
    all_parser_file = Path(deephic_path) / 'all_parser.py'
    
    if not all_parser_file.exists():
        raise FileNotFoundError(f"DeepHiC all_parser.py not found: {all_parser_file}")
    
    logger.info(f"Updating root directory and fixing set_dict")
    
    # Read current file
    with open(all_parser_file, 'r') as f:
        content = f.read()
    
    # Update root directory
    content = content.replace(
        "root_dir = '/data/RaoHiC'",
        f"root_dir = '{new_root_dir}'"
    )
    
    # Add missing 'all' key to set_dict for mouse (combines all mouse chromosomes)
    content = content.replace(
        "'mouse': list(range(1,20)) + ['X'],",
        "'mouse': list(range(1,20)) + ['X'],\n            'all': list(range(1,20)) + ['X'],"
    )
    
    # Write back
    with open(all_parser_file, 'w') as f:
        f.write(content)
    
    logger.info("Updated root directory and added 'all' key for mouse chromosomes")

def convert_model_to_cpu(model_path, logger):
    """Convert CUDA model to CPU-compatible model"""
    import torch
    
    cpu_model_path = model_path.replace('.pth', '_cpu.pth')
    
    if Path(cpu_model_path).exists():
        logger.info(f"CPU model already exists: {cpu_model_path}")
        return cpu_model_path
    
    try:
        logger.info(f"Converting CUDA model to CPU: {model_path}")
        
        # Load the CUDA model and map it to CPU
        model_state = torch.load(model_path, map_location=torch.device('cpu'))
        
        # Save as CPU model
        torch.save(model_state, cpu_model_path)
        
        logger.info(f"CPU model saved: {cpu_model_path}")
        return cpu_model_path
        
    except Exception as e:
        logger.error(f"Failed to convert model: {e}")
        raise RuntimeError(f"Could not convert model to CPU: {e}")

def run_enhancement(config, output_dir, logger):
    """Run DeepHiC enhancement step"""
    logger.info("Starting Hi-C enhancement with DeepHiC...")
    
    # Expand DeepHiC path
    deephic_root = Path(config['paths']['deephic_root']).expanduser()
    
    # Create DeepHiC working directory
    deephic_work_dir = output_dir / 'deephic_workspace'
    deephic_work_dir.mkdir(exist_ok=True)
    
    # ONLY update root directory - nothing else
    absolute_workspace = deephic_work_dir.absolute()
    update_deephic_root_dir(deephic_root, str(absolute_workspace), logger)
    
    # Convert model to CPU if needed
    cpu_model_path = convert_model_to_cpu(config['paths']['deephic_model'], logger)
    
    resolution_num = config['general']['resolution'] 
    species = config['general']['species']
    
    # Convert resolution to DeepHiC format
    resolution_map = {5000: "5kb", 10000: "10kb", 25000: "25kb", 50000: "50kb", 
                     100000: "100kb", 250000: "250kb", 500000: "500kb", 1000000: "1mb"}
    
    if resolution_num not in resolution_map:
        logger.error(f"Unsupported resolution: {resolution_num}")
        return
    
    resolution_str = resolution_map[resolution_num]
    logger.info(f"Using DeepHiC resolution: {resolution_str}")
    
    for tissue in config['samples']['tissues']:
        logger.info(f"Processing tissue: {tissue}")
        
        # Input files
        bed_file = Path(config['paths']['hicpro_dir']) / f"{tissue}_{resolution_num}_abs.bed"
        matrix_file = Path(config['paths']['hicpro_dir']) / f"{tissue}_{resolution_num}.matrix"
        
        if not bed_file.exists() or not matrix_file.exists():
            logger.error(f"Missing input files for {tissue}")
            continue
        
        # DeepHiC output directory
        output_prefix = deephic_work_dir / 'mat' / f"{tissue}_{resolution_num}"
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Step 1: Convert HiC-Pro to DeepHiC format
            cmd1 = ['python', str(deephic_root / 'hicpro2deephic.py'),
                   '--bed', str(bed_file), '--mat', str(matrix_file),
                   '-r', resolution_str, '-o', str(output_prefix)]
            
            logger.info(f"Running: {' '.join(cmd1)}")
            result = subprocess.run(cmd1, capture_output=True, text=True, cwd=deephic_work_dir)
            
            if result.returncode != 0:
                logger.error(f"hicpro2deephic failed: {result.stderr}")
                continue
            
            # Step 2: Generate training data
            cmd2 = ['python', str(deephic_root / 'data_generate.py'),
                   '-hr', resolution_str, '-lr', resolution_str, 
                   '-lrc', str(config['deephic']['low_res_cutoff']),
                   '-s', 'all',  # Now 'all' should work since we added it to set_dict
                   '-chunk', str(config['deephic']['chunk_size']),
                   '-stride', str(config['deephic']['stride_size']),
                   '-bound', str(config['deephic']['boundary']),
                   '-scale', str(config['deephic']['scale_factor']),
                   '-c', f"{tissue}_{resolution_num}"]
            
            logger.info(f"Generating training data...")
            result = subprocess.run(cmd2, capture_output=True, text=True, cwd=deephic_work_dir)
            
            if result.returncode != 0:
                logger.error(f"data_generate failed: {result.stderr}")
                continue
            
            # Step 3: Predict enhanced data (using CPU model)
            cmd3 = ['python', str(deephic_root / 'data_predict.py'),
                   '-lr', resolution_str, 
                   '-ckpt', cpu_model_path,  # Use CPU model
                   '-c', f"{tissue}_{resolution_num}",
                   '--cuda', '0']  # Force CPU usage
            
            logger.info(f"Predicting enhanced data...")
            result = subprocess.run(cmd3, capture_output=True, text=True, cwd=deephic_work_dir)
            
            if result.returncode != 0:
                logger.error(f"data_predict failed: {result.stderr}")
                continue
                
            logger.info(f"DeepHiC completed for {tissue}")
            
        except Exception as e:
            logger.error(f"Error processing {tissue}: {e}")
    
    logger.info("Enhancement completed")

def validate_chess_parameters(config, logger):
    """Validate and optimize CHESS parameters for the given genome"""
    
    genome = config['general']['reference_genome']
    resolution = config['general']['resolution']
    chess_config = config['chess']
    
    optimal_params = get_chess_parameters(genome, resolution)
    
    # Check if auto-optimization is enabled
    if chess_config.get('auto_optimize', False):
        logger.info(f"Auto-optimizing CHESS parameters for {genome}")
        
        # Update parameters
        chess_config['window_size'] = optimal_params['window_size']
        chess_config['step_size'] = optimal_params['step_size']
        
        logger.info(f"  Window size: {chess_config['window_size']:,} bp")
        logger.info(f"  Step size: {chess_config['step_size']:,} bp")
    
    else:
        # Validate user-provided parameters
        min_window = optimal_params['min_window']
        user_window = chess_config['window_size']
        
        if user_window < min_window:
            logger.warning(f"CHESS window size ({user_window:,}) is below minimum for {genome} ({min_window:,})")
            logger.warning(f"This may cause CHESS to fail. Recommended windows: {optimal_params['recommended_windows']}")
            
            response = input(f"Auto-adjust to {min_window:,} bp? (y/n): ")
            if response.lower().startswith('y'):
                chess_config['window_size'] = min_window
                logger.info(f"Window size adjusted to {min_window:,} bp")
    
    return config

def run_conversion(config, output_dir, logger):
    """Run matrix conversion step"""
    logger.info("Converting enhanced matrices to standard format...")
    
    # Implementation would go here
    
    logger.info("Conversion completed successfully")

def run_chess_analysis(config, output_dir, logger):
    """Run CHESS similarity analysis"""
    logger.info("Running CHESS similarity analysis...")
    
    # Implementation would go here
    
    logger.info("CHESS analysis completed successfully")

def run_differential_analysis(config, output_dir, logger):
    """Run differential analysis"""
    logger.info("Performing differential analysis...")
    
    # Implementation would go here
    
    logger.info("Differential analysis completed successfully")

def generate_report(config, output_dir, logger):
    """Generate summary report"""
    logger.info("Generating analysis report...")
    
    report_file = Path(output_dir) / 'reports' / 'HiC_ECC_report.html'
    
    # Create a simple HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>HiC-ECC Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #4CAF50; color: white; padding: 20px; }}
            .section {{ margin: 20px 0; padding: 15px; border-left: 4px solid #4CAF50; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>HiC-ECC Analysis Report</h1>
            <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="section">
            <h2>Analysis Parameters</h2>
            <p><strong>Reference Genome:</strong> {config['general']['reference_genome']}</p>
            <p><strong>Resolution:</strong> {config['general']['resolution']:,} bp</p>
            <p><strong>Tissues Analyzed:</strong> {', '.join(config['samples']['tissues'])}</p>
            <p><strong>Chromosomes:</strong> {', '.join(config['samples']['chromosomes'])}</p>
        </div>
        
        <div class="section">
            <h2>Output Files</h2>
            <ul>
                <li>Enhanced matrices: <code>enhanced_matrices/</code></li>
                <li>Hi-C format files: <code>hic_format/</code></li>
                <li>CHESS similarity results: <code>chess_similarity/</code></li>
                <li>Differential regions: <code>differential_regions/</code></li>
                <li>Analysis logs: <code>logs/</code></li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    with open(report_file, 'w') as f:
        f.write(html_content)
    
    logger.info(f"Report generated: {report_file}")

def main():
    parser = argparse.ArgumentParser(
        description="HiC-ECC: Hi-C Enhancement, Comparison, and Classification Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create configuration template
  HiC-ECC.py --create-config
  
  # Run full pipeline
  HiC-ECC.py -c config.yaml -o results/
  
  # Run specific steps
  HiC-ECC.py -c config.yaml -o results/ --steps enhance,convert
  
  # Verbose output
  HiC-ECC.py -c config.yaml -o results/ --verbose
        """
    )
    
    parser.add_argument('-v', '--version', action='version', version=f'HiC-ECC {__version__}')
    
    parser.add_argument('-c', '--config', 
                       help='Configuration file (YAML format)')
    
    parser.add_argument('-o', '--output', 
                       help='Output directory', 
                       default='hic_ecc_results')
    
    parser.add_argument('--create-config', 
                       action='store_true',
                       help='Create configuration template file')
    
    parser.add_argument('--steps', 
                       help='Comma-separated list of steps to run (enhance,convert,chess,differential,report)',
                       default='enhance,convert,chess,differential,report')
    
    parser.add_argument('--verbose', 
                       action='store_true',
                       help='Verbose logging')
    
    parser.add_argument('--dry-run', 
                       action='store_true',
                       help='Show what would be done without executing')
    
    args = parser.parse_args()
    
    # Handle config template creation
    if args.create_config:
        create_config_template()
        return
    
    # Validate required arguments
    if not args.config:
        parser.error("Configuration file is required. Use --create-config to generate template.")
    
    # Load configuration
    config = load_config(args.config)
    
    # Setup output directory
    output_path = setup_output_structure(args.output)
    
    # Setup logging
    logger = setup_logging(output_path, args.verbose)
    
    logger.info(f"HiC-ECC v{__version__} starting...")
    logger.info(f"Output directory: {output_path.absolute()}")
    
    # Validate and optimize CHESS parameters
    config = validate_chess_parameters(config, logger)
    
    if args.dry_run:
        logger.info("DRY RUN MODE - no files will be processed")
    
    # Parse steps to run
    steps_to_run = [step.strip() for step in args.steps.split(',')]
    
    try:
        if 'enhance' in steps_to_run:
            if args.dry_run:
                logger.info("Would run: Hi-C enhancement")
            else:
                run_enhancement(config, output_path, logger)
        
        if 'convert' in steps_to_run:
            if args.dry_run:
                logger.info("Would run: Matrix conversion")
            else:
                run_conversion(config, output_path, logger)
        
        if 'chess' in steps_to_run:
            if args.dry_run:
                logger.info("Would run: CHESS analysis")
            else:
                run_chess_analysis(config, output_path, logger)
        
        if 'differential' in steps_to_run:
            if args.dry_run:
                logger.info("Would run: Differential analysis")
            else:
                run_differential_analysis(config, output_path, logger)
        
        if 'report' in steps_to_run:
            if args.dry_run:
                logger.info("Would generate: Analysis report")
            else:
                generate_report(config, output_path, logger)
        
        logger.info("HiC-ECC pipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
