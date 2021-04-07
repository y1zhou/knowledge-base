#!/usr/bin/env python3
# %%
import re
from pathlib import Path
from random import random
from time import sleep

import requests


def featured_img_exists(content_dir: Path) -> bool:
    imgs = [x for x in content_dir.glob("featured.*")]
    return len(imgs) != 0


def read_front_matter(md_file: Path) -> str:
    if not md_file.is_file():
        raise OSError(f"File '{md_file.absolute()}' not found.'")
    with open(md_file, "r") as f:
        for _ in range(100):
            l = f.readline().strip()
            if l.startswith("unsplash_id"):
                res = re.sub(r"^unsplash_id: \"(.*)\".*$", r"\1", l)
                return res
    return ""


def retrieve_img(unsplash_id: str, dl_path: Path, ss: requests.Session) -> None:
    # sanity checks
    if not unsplash_id:
        raise ValueError(f"unsplash_id is empty.")
    if dl_path.is_file():
        raise OSError(f"{dl_path.absolute()} already exists.")

    r = ss.head(f"https://source.unsplash.com/{unsplash_id}/3000x1854")
    if r.ok:
        true_url = r.headers.get("Location")
        if true_url:
            img = ss.get(true_url)
            if img.ok:
                print(f"-- Downloading {true_url} to {dl_path}")
                with open(dl_path, "wb") as f:
                    f.write(img.content)


contents = Path("content/series")
md_files = contents.glob("**/*index.Rmarkdown")
ss = requests.session()

# %%
for md_file in md_files:
    content_dir = md_file.parent
    if featured_img_exists(content_dir):
        continue

    unsplash_id = read_front_matter(md_file)
    if unsplash_id:
        print(f"Found {unsplash_id} in {md_file}")
        img_path = content_dir / "featured.jpg"
        retrieve_img(unsplash_id, img_path, ss)
        sleep(2 * random())
