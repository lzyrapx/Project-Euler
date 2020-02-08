'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:17:23
'''
# coding: utf-8
__author__ = 'LzyRapx'
import json

def main(difficulty):
    filename = r'1_700.json'
    with open(filename) as f:
        data = json.load(f)
    # print(data)
    # print(type(data))
    cnt = 0
    print("```")
    for prob_id in data:
        # print("prob_id = ", prob_id)
        for key in data[prob_id]:
            # print("key = ", key)
            # print(data[prob_id][key])
            if key == 'info':
                if 'difficulty' in data[prob_id][key]:
                    # print("info = ", data[prob_id][key])
                    if data[prob_id][key]['difficulty'] == difficulty:
                        print(prob_id, end= " ")
                        cnt += 1
    print("\n")
    print("```")
    # print("cnt = ", cnt)
if __name__ == '__main__':
    for difficulty in range(5,101,5):
        print("difficulty = ", difficulty)
        main(difficulty)