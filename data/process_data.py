import pandas as pd

def process_single(subj_id, to_path):
    if subj_id < 10:
        path = "./subjects/0" + str(subj_id) + ".tsv"
    else:
        path = "./subjects/" + str(subj_id) + ".tsv"

    columns = ["Trial", "RiskType", "RewardType", "ResponseType", "Reward", "Shock"]

    # wiriting the header to file
    df = pd.DataFrame([], columns=columns)
    df.to_csv(to_path, mode="w+", header=True, sep="\t")

    for chunk in pd.read_csv(path, chunksize=500, sep="\t"):
        df = chunk[columns].dropna()

        df.to_csv(to_path, mode='a', header=False, sep="\t")

def process_data(to_path):

    columns = ["Trial", "RiskType", "RewardType", "ResponseType", "Reward", "Shock", "SubjID"]

    # wiriting the header to file
    df = pd.DataFrame([], columns=columns)
    df.to_csv(to_path, mode="w+", header=True, sep="\t")

    for subj_id in range(1,37):
        if subj_id == 14: continue # missing data
        if subj_id < 10:
            path = "./subjects/0" + str(subj_id) + ".tsv"
        else:
            path = "./subjects/" + str(subj_id) + ".tsv"


        for chunk in pd.read_csv(path, chunksize=500, sep="\t"):
            df = chunk[["Trial", "RiskType", "RewardType", "ResponseType", "Reward", "Shock"]].dropna()
            df["SubjID"] = subj_id
            df.to_csv(to_path, mode='a', header=False, sep="\t")

if __name__ == "__main__":
    process_data("./AllSubjectsProcessed.tsv")
