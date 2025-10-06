import easyvvuq as uq
import chaospy as cp
import matplotlib.pyplot as plt
from dask.distributed import Client
import argparse
import os

SOFTWARE_PATH='/home/nan/Documents/github/hemoflow/'
SOFTWARE_PATH=os.path.abspath(SOFTWARE_PATH)
#HEMOFLOW_PATH='/mnt/d/1_Github/VascuTreatCFD/hemoflowcfd/build/hemoflow'
HEMOFLOW_PATH=os.path.join(SOFTWARE_PATH,'build','hemoFlow')
HEMOFLOW_PATH=os.path.abspath(HEMOFLOW_PATH)
PREPRPOCESSOR_PATH=os.path.join(SOFTWARE_PATH,'preprocessor','main.py')
PREPRPOCESSOR_PATH=os.path.abspath(PREPRPOCESSOR_PATH)
#For LBMpost use the fix-vvuq branch
LBMPOST_PATH='/mnt/d/1_Github/LBMPOST/main.py'
TEMPLATE_DIR_PATH='campaign_dir'
TEMPLATE_DIR_PATH=os.path.abspath(TEMPLATE_DIR_PATH)


def run_easy_vvuq_example(client_param):
    """
        Usage of EasyVVUQ for vascular simulation.
        """
    work_dir = os.path.dirname(os.path.abspath(__file__))
    campaign = uq.Campaign(
        name="_hemoflow_sensitivity_",
        work_dir=work_dir
    )

    params = {
        "l": {"type": "float", "default": 1.0},
        "q": {"type": "float", "default": 1.5},
        "dt": {"type": "float", "default": 1e-5},
        "dx": {"type": "float", "default": 2e-4},
        "elem": {"type": "integer", "default": 20000},
        "save_dt": {"type": "integer", "default": 50},
    }

    encoder = uq.encoders.GenericEncoder(template_fname='{}/campaign_dir/input.template'.format(work_dir), delimiter='$', target_filename='input.xml')
    encoder_vox = uq.encoders.GenericEncoder(template_fname='{}/campaign_dir/vox_config.template'.format(work_dir), delimiter='$', target_filename='input/input_pipe_vox.config')
    decoder = uq.decoders.SimpleCSV(target_filename='postProc/post_data.csv', output_columns=['ane_0_velocity_vol_avg'])

    actions = uq.actions.Actions(
        uq.actions.CreateRunDirectory(root=work_dir, flatten=True),
        uq.actions.ExecuteLocal("cp --recursive "+TEMPLATE_DIR_PATH+"/. ./"),
        uq.actions.Encode(encoder_vox),
        #uq.actions.ExecuteLocal(f"conda run --live-stream -n prepost python {PREPRPOCESSOR_PATH} ./input/input_pipe_vox.config"),
        uq.actions.Encode(encoder),
        #!RUNNING ON 6 CORES
        #uq.actions.ExecuteLocal("mpirun -n 6 "+HEMOFLOW_PATH+" input.xml"), 
        #conda env for running LBMpost, livestream for stdio
        uq.actions.ExecuteLocal("conda run --live-stream -n prepost python "+LBMPOST_PATH+" ./ full"),
        uq.actions.Decode(decoder)
    )

    campaign.add_app(
        name="hemoflow_sensitivity_",
        params=params,
        actions=actions
    )

    vary = {
        "l": cp.Normal(100, 0.1),
        "q": cp.Normal(200, 0.1),
    }

    campaign.set_sampler(uq.sampling.SCSampler(vary=vary, polynomial_order=1))

    campaign.execute(pool=client_param).collate()
    print("Simulation finished.")

    print("Analysis started.")
    campaign.apply_analysis(
        uq.analysis.SCAnalysis(
            sampler=campaign.get_active_sampler(),
            qoi_cols=["ane_0_velocity_vol_avg"]
        )
    )

    results = campaign.get_last_analysis()
    plt.axis('off')

    results.plot_sobols_treemap('ane_0_velocity_vol_avg', figsize=(10, 10), filename="result_velo.png")
    import time
    time.sleep(180)


if __name__ == '__main__':

    """
    Parsing arguments to specify type of run. Possible options:
    --local - running locally.
    --slurm - running using SLURM with default option dask-jobqueue - deploy Dask 
              on common job queuing systems on HPC. Dask-mpi option also available 
              for deploying Dask from within an existing MPI environment.
    """
    parser = argparse.ArgumentParser(
        description="EasyVVUQ applied to a vertical tube deflection (using DASK).",
        epilog="",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--local", "-l", help="Run locally.", action="store_true", default=False)
    parser.add_argument("--slurm", "-s", help="Run using SLURM. Possible options: dask-jobqueue, dask-mpi.",
                        default="dask-jobqueue")

    args = parser.parse_args()

    """
    Creating client from dask.distributed according chosen run type.
    """
    if args.local:
        print("Running locally")
        client = Client(processes=True,n_workers=1, threads_per_worker=1)
        run_easy_vvuq_example(client)
    elif args.slurm == "dask-jobqueue":
        print("Running with SLURM using dask-jobqueue.")
        from dask_jobqueue import SLURMCluster

        cluster = SLURMCluster(
            shebang='#!/bin/bash -l',
            job_extra=[
                '--output /net/archive/groups/plggsano/EasyVVUQ/slurm_outputs/slurm-%j.out',
                '--error /net/archive/groups/plggsano/EasyVVUQ/slurm_outputs/slurm-%j.err'],
            queue='plgrid-short',
            project='plgsano2',
            cores=24,
            processes=8,
            memory="20GB",
            walltime='00:10:00',
            interface='ib0',
            scheduler_options={'interface': 'eth55'}
        )
        cluster.job_cls.submit_command = 'sbatch'
        cluster.submit_command = 'sbatch'
        cluster.scale(jobs=8)
        client = Client(cluster)
        print(cluster)
        print(client)
        run_easy_vvuq_example(client)
    elif args.slurm == "dask-mpi":
        print("Running with SLURM using dask-mpi.")
        from dask_mpi import initialize
        from distributed.scheduler import logger
        import socket

        initialize()
        client = Client()

        host = client.run_on_scheduler(socket.gethostname)
        port = client.scheduler_info()['services']['dashboard']
        login_node_address = "login@pro.cyfronet.pl"  # Change this to the address/domain of your login node

        logger.info(f"ssh -N -L {port}:{host}:{port} {login_node_address}")

        from dask.distributed import performance_report

        with performance_report(filename="dask-report.html"):
            run_easy_vvuq_example(client)

        client.profile(filename="dask-profile.html")

    else:
        print("Incorrect slurm option specified!")
        exit(1)
