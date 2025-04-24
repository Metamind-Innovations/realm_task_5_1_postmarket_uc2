import kfp
from kfp import dsl
from kfp.dsl import Output, Model, Dataset


@kfp.dsl.component(
    base_image="python:3.10",
    packages_to_install=["requests"]
)
def download_project(
        github_repo_url: str,
        project_files: Output[Model],
        real_world_data: Output[Dataset],
        synthetic_data: Output[Dataset],
        groundtruth_labels: Output[Dataset],
):
    import subprocess
    import shutil
    import re
    from pathlib import Path

    branch_match = re.search(r"/tree/([^/]+)", github_repo_url)
    if branch_match:
        url_branch = branch_match.group(1)
        print(f"Extracted branch '{url_branch}' from URL")
        branch = url_branch
    repo_url = re.sub(r"/tree/[^/]+/?$", "", github_repo_url.strip())
    if not repo_url.endswith(".git"):
        repo_url = repo_url.rstrip("/") + ".git"
    print(f"Original URL: {github_repo_url}")
    print(f"Using branch: {branch}")
    print(f"Using repository URL: {repo_url}")

    print("Installing git...")
    subprocess.run(["apt-get", "update"], check=True, capture_output=True, text=True)
    subprocess.run(["apt-get", "install", "-y", "git"], check=True, capture_output=True, text=True)
    print("Git installed.")
    temp_dir = Path("/tmp/pharmcat_repo")
    if temp_dir.exists():
        print(f"Removing existing temp directory: {temp_dir}")
        shutil.rmtree(temp_dir)
    clone_cmd = ["git", "clone", "--depth", "1", "-b", branch, repo_url, str(temp_dir)]
    print(f"Running: {' '.join(clone_cmd)}")
    try:
        subprocess.run(clone_cmd, check=True, capture_output=True, text=True)
        print(f"Successfully cloned branch '{branch}' to {temp_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Git clone failed with error: {e.stderr}")
        raise Exception(f"Failed to clone repository: {repo_url}, branch: {branch}")

    project_files_path = Path(project_files.path)
    real_world_data_path = Path(real_world_data.path)
    synthetic_data_path = Path(synthetic_data.path)
    groundtruth_labels_path = Path(groundtruth_labels.path)

    project_files_path.mkdir(parents=True, exist_ok=True)
    real_world_data_path.mkdir(parents=True, exist_ok=True)
    synthetic_data_path.mkdir(parents=True, exist_ok=True)
    groundtruth_labels_path.mkdir(parents=True, exist_ok=True)
    print(f"Created KFP output directories.")

    items_to_copy = {
        "project_files": [
            "statistical_analysis.py",
            "expert_knowledge.py",
            "adversarial_evaluation.py",
            "requirements.txt",
        ],
        "real_world_data": ["data/rwd_vcf/"],
        "synthetic_data": ["data/synthetic_vcf/"],
        "groundtruth_labels": ["data/groundtruth/"],
    }
    found_items = {key: [] for key in items_to_copy}
    copied_items_count = 0

    print(f"\nSearching for required items in cloned repo: {temp_dir}")
    for target_list_key, target_items in items_to_copy.items():
        for item_name in target_items:
            src_path = temp_dir / item_name
            dst_path = None
            is_dir = src_path.is_dir()
            is_file = src_path.is_file()

            if not (is_dir or is_file):
                print(f"  Warning: Item '{item_name}' not found at {src_path}")
                continue

            try:
                if target_list_key == "project_files":
                    dst_path = project_files_path / item_name
                    if is_dir:
                        print(f"Copying directory '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                    else:
                        print(f"Copying file '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copy2(src_path, dst_path)

                elif target_list_key == "real_world_data":
                    dst_path = real_world_data_path
                    if is_dir:
                        print(f"Copying directory '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                    else:
                        print(f"Copying file '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copy2(src_path, dst_path)

                elif target_list_key == "synthetic_data":
                    dst_path = synthetic_data_path
                    if is_dir:
                        print(f"Copying directory '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                    else:
                        print(f"Copying file '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copy2(src_path, dst_path)

                elif target_list_key == "groundtruth_labels":
                    dst_path = groundtruth_labels_path
                    if is_dir:
                        print(f"Copying directory '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                    else:
                        print(f"Copying file '{item_name}' from {src_path} to {dst_path}...")
                        shutil.copy2(src_path, dst_path)

                if dst_path:
                    found_items[target_list_key].append(item_name)
                    copied_items_count += 1
                    print(f"  Successfully copied '{item_name}' to {dst_path}")

            except Exception as e:
                print(f"  Error copying {item_name} from {src_path} to {dst_path}: {e}")

    print(f"\nFinished search. Copied {copied_items_count} items.")

    print(f"Removing temporary clone directory: {temp_dir}")
    shutil.rmtree(temp_dir)
    print("Download component finished.")

    return (project_files, real_world_data, synthetic_data, groundtruth_labels)


@dsl.component(base_image="python:3.10")
def expert_knowledge(
    project_files: dsl.Input[dsl.Model],
    input_dir: dsl.Input[dsl.Dataset],
    output_file: dsl.Output[dsl.Artifact],
):
    """Execute expert_knowledge.py with specific arguments"""
    import subprocess
    import os
    from pathlib import Path

    script_path = Path(project_files.path) / "expert_knowledge.py"
    
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found at {script_path}")

    command = [
        "python",
        str(script_path),
        "--input_dir",
        input_dir.path,
        "--output_file",
        output_file.path,
    ]

    print(f"\nRunning expert knowledge evaluation: {' '.join(command)}")

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("Expert knowledge evaluation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during expert knowledge evaluation: {e.stdout}\n{e.stderr}")
        raise Exception("Expert knowledge evaluation failed.")


@dsl.component(base_image="python:3.10")
def statistical_analysis(
    project_files: dsl.Input[dsl.Model],
    input_dir: dsl.Input[dsl.Dataset],
    output_file: dsl.Output[dsl.Artifact],
):
    """Execute statistical_analysis.py with specific arguments"""
    import subprocess
    from pathlib import Path

    script_path = Path(project_files.path) / "statistical_analysis.py"
    
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found at {script_path}")

    command = [
        "python",
        str(script_path),
        "--input_dir",
        input_dir.path,
        "--output_file",
        output_file.path,
    ]

    print(f"\nRunning statistical analysis: {' '.join(command)}")

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("Statistical analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during statistical analysis: {e.stdout}\n{e.stderr}")
        raise Exception("Statistical analysis failed.")


@dsl.container_component
def pharmcat_analysis_docker(
    input_folder: dsl.Input[dsl.Dataset],
    result_folder: dsl.Output[dsl.Dataset],
):
    command_str = f"mkdir -p '{result_folder.path}' && python3 -u /scripts/pharmcat_folder_processor.py --input_folder '{input_folder}' --result_folder '{result_folder.path}'"
    return dsl.ContainerSpec(
        image="gigakos/pharmcat-realm:latest",  # Insert your Docker image here (e.g. "docker.io/<username>/pharmcat-realm:latest")
        command=["sh", "-c"],
        args=[command_str],
    )


@dsl.component(base_image="python:3.10")
def adversarial_evaluation(
    project_files: dsl.Input[dsl.Model],
    groundtruth_file: dsl.Input[dsl.Dataset],
    rwd_predictions_file: dsl.Input[dsl.Dataset],
    synthetic_predictions_file: dsl.Input[dsl.Dataset],
    output_file: dsl.Output[dsl.Dataset],
):
    from pathlib import Path
    import subprocess

    script_path = Path(project_files.path) / "adversarial_evaluation.py"
    
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found at {script_path}")

    command = [
        "python",
        str(script_path),
        "--groundtruth_file",
        groundtruth_file.path,
        "--rwd_predictions_file",
        rwd_predictions_file.path,
        "--synthetic_predictions_file",
        synthetic_predictions_file.path,
        "--output_file",
        output_file.path,
    ]

    print(f"\nRunning adversarial evaluation: {' '.join(command)}")

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("Adversarial evaluation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during adversarial evaluation: {e.stdout}\n{e.stderr}")
        raise Exception("Adversarial evaluation failed.")


@dsl.pipeline(
    name="Post-Market Evaluation Pipeline",
    description="Pipeline that runs the full post-market evaluation",
)
def post_market_evaluation_pipeline(
    github_repo_url: str,
):
    
    # Step 0: Download the project
    download_task = download_project(
        github_repo_url=github_repo_url,
    )
    download_task.set_caching_options(False)

    # Step 1: Expert Knowledge
    expert_task = expert_knowledge(
        project_files=download_task.outputs["project_files"],
        input_dir=download_task.outputs["synthetic_data"],
    )

    # Step 2: Statistical Analysis
    stats_task = statistical_analysis(
        project_files=download_task.outputs["project_files"],
        input_dir=download_task.outputs["synthetic_data"], 
    )

    # Step 3a: PharmCAT Analysis for RWD
    pharmcat_task_rwd = pharmcat_analysis_docker(
        input_folder=download_task.outputs["real_world_data"],
    )

    # Step 3b: PharmCAT Analysis for Synthetic Data
    pharmcat_task_synthetic = pharmcat_analysis_docker(
        input_folder=download_task.outputs["synthetic_data"],
    )

    # Step 4: Adversarial Evaluation
    adv_task = adversarial_evaluation(
        project_files=download_task.outputs["project_files"],
        groundtruth_file=download_task.outputs["groundtruth_labels"],
        rwd_predictions_file=pharmcat_task_rwd.outputs["result_folder"],
        synthetic_predictions_file=pharmcat_task_synthetic.outputs["result_folder"],
    )


if __name__ == "__main__":
    kfp.compiler.Compiler().compile(
        pipeline_func=post_market_evaluation_pipeline,
        package_path="post_market_evaluation_pipeline.yaml",
    )
