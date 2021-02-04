import json
import os

def parse_count ():
    path = "/home/jamie/nthu_cs/research/projects/gci/project/output/construct"
    for entry in os.scandir(path):
        if entry.is_file():
            os.rename(entry.path, os.path.splitext(entry.path)[0] + ".gci.json")
            # with open(entry.path) as json_file:
            #     json_tree = json.load(json_file)
            #     print(json_tree["timeEnd"]-json_tree["timeStart"])

if __name__ == "__main__":
    parse_count()
