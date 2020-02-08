'''
@Author: zhaoyang.liang
@Github: https://github.com/LzyRapx
@Date: 2020-02-01 16:17:23
'''
#coding: utf-8
__author__ = 'LzyRapx'

import requests
import re
import json
import shutil
from bs4 import BeautifulSoup

def get_title(soup):
    return soup.h2.text

def get_info(soup):
    # get raw info section
    raw_info = soup.find('span', {'style': 'left:-250px;width:250px;top:30px;'})
    # print("raw info = ", raw_info)
    # get published, solved, and difficulty information using regex
    published = re.search(r'Published on (.+);', str(raw_info))
    solved = re.search(r'Solved by (\d+)', str(raw_info))
    difficulty = re.search(r'Difficulty rating: (\d+)%', str(raw_info))
    return {
        'published': published.group(1),
        'solved': int(solved.group(1)),
        # since newly published problems don't have a difficulty, put None if it isn't there
        'difficulty': int(difficulty.group(1)) if difficulty != None else None
    }
def get_images(problem_content):
    images = []
    for tag in problem_content.descendants:
        if tag.name == 'img' and tag['src'] != 'images/spacer.gif':
            path = tag['src']
            images.append(path)
            r = requests.get('http://projecteuler.net/' + path, stream=True)
            with open('./data/images/{}'.format(path[15:]), 'wb+') as f:
                shutil.copyfileobj(r.raw, f)
    return images

def get_files(problem_content):
    files = []
    for tag in problem_content.descendants:
        if tag.name == 'a' and 'project' in tag['href']:
            path = tag['href']
            files.append(path)

            r = requests.get('http://projecteuler.net/' + path)
            with open('./data/files/{}'.format(path[18:]), 'w+') as f:
                f.write(r.text)
    return files

def get_content(problem_content):
    # somehow unwrapping the <div> tag would be better, can't get it to work though
    # so I'll just use regex to get stuff inside the tag
    html = re.search(r'<div class="problem_content" role="problem">(.*)<\/div>', str(problem_content), re.S).group(1)
    # get just text
    text = problem_content.get_text()
    return {
        'text': text,
        'html': html,
        'images': get_images(problem_content),
        'files': get_files(problem_content)
    }

def scrape(problem):
    url = 'https://projecteuler.net/problem=' + str(problem)
    r = requests.get(url)
    raw_html = r.text
    soup = BeautifulSoup(raw_html, 'html.parser')
    # get raw problem content
    problem_content = soup.find('div', {'class': 'problem_content'})

    return {
        'number': problem,
        'url': url,
        'title': get_title(soup),
        'info': get_info(soup),
        'content': get_content(problem_content)
    }

def get_num_problems():
    r = requests.get('http://projecteuler.net/recent')
    raw_html = r.text
    soup = BeautifulSoup(raw_html, 'html.parser')
    table_rows = soup.find('table', {'id': 'problems_table'}).find_all('tr')
    most_recent = table_rows[1]
    most_recent_number = most_recent.find('td').text
    return int(most_recent_number)

if __name__ == '__main__':
    start = 1
    stop = get_num_problems()
    print("total = ", stop)
    filename = './data/{}_{}.json'.format(start, stop)
    output = {}
    print("start")
    for num in range(start, stop + 1):
        print('Scraping problem {} / {}'.format(num, stop))
        output[num] = scrape(num)
    with open(filename, 'w+') as f:
        f.write(json.dumps(output, indent=2))
    print('\nProblems {} - {} successfully scraped and saved to {}'.format(start, stop, filename))