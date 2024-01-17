import os
import argparse
import requests
from tqdm import tqdm
import shutil
from concurrent.futures import ThreadPoolExecutor
import subprocess


def download_file(url, filename, position):
    """使用requests下载文件并用tqdm显示进度条"""
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    chunk_size = 1024
    with open(filename, 'wb') as f, tqdm(
        desc=filename,
        total=total_size,
        unit='iB',
        unit_scale=True,
        position=position,
        leave=False  # 设置为False，下载完成后隐藏进度条
    ) as bar:
        for chunk in response.iter_content(chunk_size=chunk_size):
            size = f.write(chunk)
            if size is not None:
                bar.update(int(size))

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='dbCAN Database Downloader')
    parser.add_argument('--cpus', type=int, default=8, help='Number of CPUs for parallel downloads')
    parser.add_argument('--db-dir', type=str, default='db', help='Database directory name')
    parser.add_argument('--clean', action='store_true', help='Clean up temporary files after processing')
    return parser.parse_args()
def clean_up(db_dir):
    try:
        shutil.rmtree(db_dir)
        print(f"Folder '{db_dir}' and all its contents have been deleted.")
    except OSError as error:
        print(f"Error: {error}")
        print(f"Failed to delete the folder '{db_dir}'.")
def get_remote_file_size(url):
    """获取远程文件的大小"""
    response = requests.head(url, allow_redirects=True)
    if 'content-length' in response.headers:
        return int(response.headers['content-length'])
    else:
        return 0
def run_command(cmd):
    """运行单个shell命令"""
    subprocess.run(cmd, shell=True, check=True)

def main():
    args = parse_args()

    if args.clean and os.path.exists(args.db_dir):
        clean_up(args.db_dir)

    if not os.path.exists(args.db_dir):
        os.makedirs(args.db_dir)
    os.chdir(args.db_dir)

    download_urls = [
    "https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv",
    "https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa",
    "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx",
    "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz",
    "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm",
    "https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm"
    ]

    total_size = sum(get_remote_file_size(url) for url in download_urls)

    with ThreadPoolExecutor(max_workers=args.cpus) as executor:
        futures = []
        for i, url in enumerate(download_urls):
            filename = os.path.basename(url)
            futures.append(executor.submit(download_file, url, filename, i))
        
        # 显示总的进度条
        with tqdm(total=total_size, unit='iB', unit_scale=True, position=len(download_urls), leave=False) as bar:
            for future in futures:
                result = future.result()
                if result is not None:
                    bar.update(int(result))

    # 执行其他命令
    other_commands = [
        "mv fam-substrate-mapping-08012023.tsv fam-substrate-mapping.tsv",
        "makeblastdb -in PUL.faa -dbtype prot",
        "mv dbCAN-PUL_12-12-2023.xlsx dbCAN-PUL.xlsx",
        "tar xzf dbCAN-PUL.tar.gz",
        "hmmpress -f dbCAN_sub.hmm",
        "mv CAZyDB.07262023.fa CAZyDB.fa",
        "diamond makedb --in CAZyDB.fa -d CAZy",
        "mv dbCAN-HMMdb-V12.txt dbCAN.txt",
        "hmmpress dbCAN.txt",
        "diamond makedb --in tcdb.fa -d tcdb",
        "hmmpress -f tf-1.hmm",
        "hmmpress -f tf-2.hmm",
        "hmmpress -f stp.hmm"
    ]
    for cmd in other_commands:
        subprocess.run(cmd, shell=True)
    # for cmd in other_commands:
    #     os.system(cmd)

    os.chdir('..')
    # 下载样本文件
    sample_commands = [
        "wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna",
        "wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa",
        "wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff"
    ]

    with ThreadPoolExecutor() as executor:
        executor.map(run_command, sample_commands)

if __name__ == "__main__":
    main()
