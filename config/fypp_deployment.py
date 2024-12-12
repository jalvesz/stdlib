import os
import fypp
import argparse
from joblib import Parallel, delayed

C_PREPROCESSED = (
    "stdlib_linalg_constants",
    "stdlib_blas",
    "stdlib_linalg_blas_aux",
    "stdlib_blas_level1",
    "stdlib_blas_level2_gen",
    "stdlib_blas_level2_ban",
    "stdlib_blas_level2_pac",
    "stdlib_blas_level2_sym",
    "stdlib_blas_level2_tri",
    "stdlib_blas_level3_gen",
    "stdlib_blas_level3_sym",
    "stdlib_blas_level3_tri",

    "stdlib_lapack",
    "stdlib_linalg_lapack_aux",
    "stdlib_lapack_auxiliary_others",
    "stdlib_lapack_auxiliary_parameters",
    "stdlib_lapack_blas_like_base",
    "stdlib_lapack_blas_like_l1",
    "stdlib_lapack_blas_like_l2",
    "stdlib_lapack_blas_like_l3",
    "stdlib_lapack_blas_like_mnorm",
    "stdlib_lapack_blas_like_scalar",
    "stdlib_lapack_cosine_sine",
    "stdlib_lapack_cosine_sine2",
    "stdlib_lapack_eigv_comp",
    "stdlib_lapack_eigv_comp2",
    "stdlib_lapack_eigv_gen",
    "stdlib_lapack_eigv_gen2",
    "stdlib_lapack_eigv_gen3",
    "stdlib_lapack_eigv_std_driver",
    "stdlib_lapack_eigv_svd_bidiag_dc",
    "stdlib_lapack_eigv_svd_drivers",
    "stdlib_lapack_eigv_svd_drivers2",
    "stdlib_lapack_eigv_svd_drivers3",
    "stdlib_lapack_eigv_sym_comp",
    "stdlib_lapack_eigv_sym",
    "stdlib_lapack_eigv_tridiag",
    "stdlib_lapack_eigv_tridiag2",
    "stdlib_lapack_eigv_tridiag3",
    "stdlib_lapack_givens_jacobi_rot",
    "stdlib_lapack_householder_reflectors",
    "stdlib_lapack_lsq",
    "stdlib_lapack_lsq_aux",
    "stdlib_lapack_lsq_constrained",
    "stdlib_lapack_orthogonal_factors_ql",
    "stdlib_lapack_orthogonal_factors_qr",
    "stdlib_lapack_orthogonal_factors_rz",
    "stdlib_lapack_others",
    "stdlib_lapack_solve_aux",
    "stdlib_lapack_solve_chol_comp",
    "stdlib_lapack_solve_chol",
    "stdlib_lapack_solve_ldl_comp",
    "stdlib_lapack_solve_ldl_comp2",
    "stdlib_lapack_solve_ldl_comp3",
    "stdlib_lapack_solve_ldl_comp4",
    "stdlib_lapack_solve_ldl",
    "stdlib_lapack_solve_lu_comp",
    "stdlib_lapack_solve_lu",
    "stdlib_lapack_solve_tri_comp",
    "stdlib_lapack_svd_bidiag_qr",
    "stdlib_lapack_svd_comp",
    "stdlib_lapack_svd_comp2",
)

def pre_process_fypp(args):
    """use fypp to preprocess all source files. 

    Processed files will be dumped at <current_folder>/temp/<file.f90> or <file.F90>

    Parameters
    ----------
    args : 
        CLI arguments.
    """
    kwd = []
    kwd.append("-DMAXRANK="+str(args.maxrank))
    kwd.append("-DPROJECT_VERSION_MAJOR="+str(args.vmajor))
    kwd.append("-DPROJECT_VERSION_MINOR="+str(args.vminor))
    kwd.append("-DPROJECT_VERSION_PATCH="+str(args.vpatch))
    if args.with_qp:
        kwd.append("-DWITH_QP=True")
    if args.with_xdp:
        kwd.append("-DWITH_XDP=True")

    optparser = fypp.get_option_parser()
    options, leftover = optparser.parse_args(args=kwd)
    options.includes = ['include']
    if args.lnumbering:
        options.line_numbering = True
    tool = fypp.Fypp(options)

    # Check destination folder for preprocessing. if not 'stdlib-fpm', it is assumed to be the root folder.
    if not os.path.exists('src'+os.sep+'temp'):
        os.makedirs('src'+os.sep+'temp')
    if not os.path.exists('test'+os.sep+'temp'):
        os.makedirs('test'+os.sep+'temp')

    # Define the folders to search for *.fypp files
    folders = ['src','test']
    # Process all folders
    fypp_files = [os.path.join(root, file) for folder in folders
                  for root, _, files in os.walk(folder)
                  for file in files if file.endswith(".fypp")]
    
    def process_f(file):
        source_file = file
        root = os.path.dirname(file)
        if not os.path.exists(root+os.sep+'temp'):
            os.makedirs(root+os.sep+'temp')
        basename = os.path.splitext(os.path.basename(source_file))[0]
        sfx = 'f90' if basename not in C_PREPROCESSED else 'F90'
        target_file = root+os.sep+'temp' + os.sep + basename + '.' + sfx
        tool.process_file(source_file, target_file)
    
    Parallel(n_jobs=args.njob)(delayed(process_f)(f) for f in fypp_files)
    
    return


def deploy_stdlib_fpm():
    """create the stdlib-fpm folder for backwards compatibility (to be deprecated)
    """
    import shutil
    prune=(
        "test_hash_functions.f90",
        "f18estop.f90",
    )
    if not os.path.exists('stdlib-fpm'+os.sep+'src'):
        os.makedirs('stdlib-fpm'+os.sep+'src')
    if not os.path.exists('stdlib-fpm'+os.sep+'test'):
        os.makedirs('stdlib-fpm'+os.sep+'test')
    if not os.path.exists('stdlib-fpm'+os.sep+'example'):
        os.makedirs('stdlib-fpm'+os.sep+'example')

    def recursive_copy(folder):
        for root, _, files in os.walk(folder):
            for file in files:
                if file not in prune:
                    if file.endswith(".f90") or file.endswith(".F90") or file.endswith(".dat") or file.endswith(".npy"):
                        shutil.copy2(os.path.join(root, file), 'stdlib-fpm'+os.sep+folder+os.sep+file)
    recursive_copy('src')
    recursive_copy('test')
    recursive_copy('example')
    for file in ['.gitignore','fpm.toml','LICENSE','VERSION']:
        shutil.copy2(file, 'stdlib-fpm'+os.sep+file)
    return

def fpm_build(args,unknown):
    import subprocess
    #==========================================
    # check compilers
    FPM_FC  = os.environ['FPM_FC']  if "FPM_FC"  in os.environ else "gfortran"
    FPM_CC  = os.environ['FPM_CC']  if "FPM_CC"  in os.environ else "gcc"
    FPM_CXX = os.environ['FPM_CXX'] if "FPM_CXX" in os.environ else "gcc"
    #==========================================
    # Filter out flags
    preprocessor = { 'gfortran':'-cpp ' , 'ifort':'-fpp ' , 'ifx':'-fpp ' }
    flags = preprocessor[FPM_FC]
    for idx, arg in enumerate(unknown):
        if arg.startswith("--flag"):
            flags= flags + unknown[idx+1]
    #==========================================
    # build with fpm
    subprocess.run("fpm build"+
                   " --compiler "+FPM_FC+
                   " --c-compiler "+FPM_CC+
                   " --cxx-compiler "+FPM_CXX+
                   " --flag \"{}\"".format(flags), shell=True, check=True)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocess stdlib source files.')
    # fypp arguments
    parser.add_argument("--vmajor", type=int, default=0, help="Project Version Major")
    parser.add_argument("--vminor", type=int, default=5, help="Project Version Minor")
    parser.add_argument("--vpatch", type=int, default=0, help="Project Version Patch")

    parser.add_argument("--njob", type=int, default=4, help="Number of parallel jobs for preprocessing")
    parser.add_argument("--maxrank",type=int, default=4, help="Set the maximum allowed rank for arrays")
    parser.add_argument("--with_qp",action='store_true', help="Include WITH_QP in the command")
    parser.add_argument("--with_xdp",action='store_true', help="Include WITH_XDP in the command")
    parser.add_argument("--lnumbering",action='store_true', help="Add line numbering in preprocessed files")
    parser.add_argument("--deploy_stdlib_fpm",action='store_true', help="create the stdlib-fpm folder")
    # external libraries arguments
    parser.add_argument("--build", action='store_true', help="Build the project")

    args, unknown = parser.parse_known_args()
    #==========================================
    # read current manifest
    with open('VERSION', 'r') as file:
        version = file.read().split(".")
        vmajor, vminor, vpatch = [int(value) for value in version]
        args.vmajor = max(vmajor,args.vmajor)
        args.vminor = max(vminor,args.vminor)
        args.vpatch = max(vpatch,args.vpatch)
    #==========================================
    # pre process the meta programming fypp files
    pre_process_fypp(args)
    if args.deploy_stdlib_fpm:
        deploy_stdlib_fpm()
    #==========================================
    # build using fpm
    if args.build:
        fpm_build(args,unknown)