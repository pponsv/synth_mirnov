import h5py
import ast


def get_keys(filename):
    keys = []
    with h5py.File(filename, "r") as hfile:
        for key in hfile:
            if key != "info":
                keys.append(key)
    return keys


def get_info(filename: str):
    with h5py.File(filename, "r") as hfile:
        time = hfile["info"]["time"][:]  # type: ignore
        name = hfile["info"]["db_name"][()].decode("utf-8")  # type: ignore
        len_th = hfile["info"]["LEN_TH"][()]  # type: ignore
        len_ph = hfile["info"]["LEN_PH"][()]  # type: ignore
    return time, name, len_th, len_ph


def get_data(filename: str, key: str):
    with h5py.File(filename, "r") as hfile:
        m = hfile[key]["m"][()]  # type: ignore
        n = hfile[key]["n"][()]  # type: ignore
        f0 = hfile[key]["f0"][()]  # type: ignore
        s0 = hfile[key]["s0"][()]  # type: ignore
        sigma = hfile[key]["sigma"][()]  # type: ignore
        alt_name = hfile[key]["alt_name"][()].decode("utf-8")  # type: ignore
        db = hfile[key]["data"][:]  # type: ignore
    return (m, n, f0, s0, sigma, alt_name), db


def get_data_new(filename: str, key: str):
    with h5py.File(filename, "r") as hfile:
        run_params = ast.literal_eval(hfile[key]["run_params"][()].decode("utf-8"))  # type: ignore
        db = hfile[key]["data"][:]  # type: ignore
    return run_params, db


def write_info(filename, time, db_name, len_th, len_ph):
    with h5py.File(filename, "w") as hfile:  #   Deletes the file if exists
        group = hfile.create_group("info")
        group.create_dataset("time", data=time)
        group.create_dataset("db_name", data=f"{len_th}_{len_ph}_{db_name}")
        group.create_dataset("LEN_TH", data=len_th)
        group.create_dataset("LEN_PH", data=len_ph)


def write_to_h5(filename, db, run_params):
    with h5py.File(filename, "a") as hfile:
        m, n, f0, s0, sigma, alt_name = run_params
        gname = f"m_{m}_n_{n}_s0_{s0:.3f}_sigma_{sigma:.4f}_{alt_name}"
        group = hfile.create_group(gname)
        group.create_dataset("n", data=n)
        group.create_dataset("m", data=m)
        group.create_dataset("f0", data=f0)
        group.create_dataset("s0", data=s0)
        group.create_dataset("sigma", data=sigma)
        group.create_dataset("alt_name", data=alt_name)
        group.create_dataset("data", data=db)


def write_to_h5_new(filename, db, key, run_params):
    with h5py.File(filename, "a") as hfile:
        group = hfile.create_group(key)
        group.create_dataset("run_params", data=str(run_params))
        group.create_dataset("data", data=db)
