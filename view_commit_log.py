#!/usr/bin/env python

program_description="Run with option -h to view help. Script to generate a report of recent commits. If the appropriate LaTeX packages are installed, can also be used to generate slides using beamer. Requires python-git, python-gitdb, python-bs4, python-feedparser and (probably) Python 2.7. Suggested LaTeX packages: latex-beamer texlive"

__author__="Kyle Barlow"

import argparse
import subprocess
from subprocess import PIPE
from git import *
import gitdb
import os
import feedparser
from bs4 import BeautifulSoup
import urllib
import cPickle as pickle
import time
import sys
from datetime import datetime
import codecs

test_server_rss='http://test.rosettacommons.org/rss'
first_git_rev_id=55120 # Should not need to be changed (ever?)
cache_pickle='.view_commit_log.pickle'
github_compare='https://github.com/RosettaCommons/main/compare/'
test_server_url='http://test.rosettacommons.org/rev/'

class Reporter:
    def __init__(self,task,report_interval=1):
        self.report_interval=report_interval # (seconds)
        self.start=time.time()
        self.lastreport=self.start
        self.task=task
        print 'Starting '+task
    def report(self,n):
        t=time.time()
        if self.lastreport<(t-self.report_interval):
            self.lastreport=t
            sys.stdout.write("  Processed: "+str(n)+" \r" )
            sys.stdout.flush()
    def done(self):
        print 'Done %s, took %.3f seconds\n' % (self.task,time.time()-self.start)

class RosettaRevision:
    def __init__(self,rev_id):
        self.rev_id=rev_id
        self.get_sha1() # Fetch corresponding hash from testing server

        self.author=None
        self.status=None
        
        self.interesting=None
        self.commit_message=None
        self.parents=None
        self.short_sha1=None
        self.date_time=None

    def __repr__(self):
        return '%d->%s'%(self.rev_id,self.sha1)

    def get_info_from_repo(self,repo):
        # Check if commit is in repo
        try:
            commit=repo.commit(self.sha1)
        except gitdb.exc.BadObject:
            print "Couldn't find sha1 %s, trying fetching" % self.sha1
            repo.remotes.origin.fetch()

        try:
            commit=repo.commit(self.sha1)
        except gitdb.exc.BadObject:
            print "Failed. Try fetching from origin and attempting again..."
            raise Exception('Try fetching')

        self.get_commit_message(commit)
        self.get_parents(commit)
        self.get_author(commit)
        self.get_date_time(commit)
        self.get_short_sha1(repo)

    def get_short_sha1(self,repo):
        self.short_sha1=repo.git.rev_parse('--short',self.sha1)

    def get_date_time(self,commit):
        self.date_time=datetime.fromtimestamp(commit.committed_date)

    def get_author(self,commit):
        self.author=str(commit.author)

    def get_parents(self,commit):
        self.parents=[commit.hexsha for commit in commit.parents]

    def get_commit_message(self,commit):
        self.commit_message=commit.message

    def get_sha1(self):
        self.sha1=None

        if self.rev_id < first_git_rev_id:
            raise Exception('Revision id %d corresponds to an SVN commit'%(self.rev_id))

        soup = BeautifulSoup(urllib.urlopen(test_server_url+'%d'%(self.rev_id)).read())
        links = soup.find_all('a')

        if not len(links) > 0:
            raise Exception('Test server page not parsed correctly')

        for link in links:
            href = link.attrs['href']
            if 'github.com/RosettaCommons/main/commit' in href:
                self.sha1 = get_hash_from_github_url(href)

        if self.sha1==None:
            raise Exception("Couldn't look up hash for revision id: %d"%(self.rev_id))

def get_hash_from_github_url(github_url):
    split_url=github_url.split('/')
    return split_url[-1]

def input_yes_no(msg=''):
    """
    Simple helper function
    """
    print '\n'+msg
    while(True):
        i=raw_input('Input yes or no: ')
        i=i.lower()
        if i=='y' or i=='yes':
            return True
        elif i=='n' or i=='no':
            return False
        else:
            print 'ERROR: Bad input. Must enter y/n/yes/no'

def check_negative(value):
    ivalue = int(value)
    if ivalue < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def check_date(value):
    try:
        dvalue = datetime.strptime(value,'%Y-%m-%d')
    except Exception:
         raise argparse.ArgumentTypeError("%s is not a date in format YYYY-MM-DD" % value)
    return dvalue

def main():

    description = program_description

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-m','--mark_interesting',
                        default=False,
                        action='store_true',
                        help='Set this flag to go through all unreviewed commits to mark which ones are interesting')

    parser.add_argument('-a','--begin_rev_id',
                        default=first_git_rev_id+1, # So that parent can be looked up
                        type=check_negative,
                        help='Test server revision ID of oldest commit to output')
    parser.add_argument('-d','--begin_date',
                        default=None,
                        type=check_date,
                        help='A date in format YYYY-MM-DD. If this option is specified, only commits on or after this date will be output')

    parser.add_argument('-e','--end_date',
                        default=None,
                        type=check_date,
                        help='A date in format YYYY-MM-DD. Only commits before (not on) this date will be output')
    parser.add_argument('-b','--end_rev_id',
                        default=None,
                        type=check_negative,
                        help='Test server revision ID of newest commit to output')

    parser.add_argument('-o','--output_slides_pdf',
                        default=None,
                        help='Location to output slides using beamer/latex (.pdf will be appended automatically)')

    parser.add_argument('-i','--interesting_only',
                        default=False,
                        action='store_true',
                        help='Only output commits marked as "interesting"')

    parser.add_argument('-g','--git_repository',
                        default='../main',
                        help='Location of git repository. Defaults to ../main')

    args=parser.parse_args()

    # Check input git directory
    assert(os.path.isdir(args.git_repository))

    # Load pickle or initialize
    if os.path.isfile(cache_pickle):
        with open(cache_pickle,'r') as f:
            rev_dict=pickle.load(f)
    else:
        rev_dict=initialize_data()

    rev_dict=update_cache(rev_dict,args.git_repository)
    
    # Program flow control
    if args.mark_interesting:
        mark_interesting(rev_dict)
    
    if args.output_slides_pdf:
        output_revisions=[]
        for revision in rev_dict:
            rev=rev_dict[revision]
            output_rev=True

            if args.interesting_only:
                if not rev.interesting:
                    output_rev=False

            if args.begin_rev_id:
                if rev.rev_id<args.begin_rev_id:
                    output_rev=False

            if args.end_rev_id:
                if rev.rev_id>args.end_rev_id:
                    output_rev=False

            if args.begin_date:
                if rev.date_time<args.begin_date:
                    output_rev=False

            if args.end_date:
                if rev.date_time>args.end_date:
                    output_rev=False

            if output_rev:
                output_revisions.append(revision)

        output_slides(args.output_slides_pdf,rev_dict,output_revisions)
        
    # Save pickle
    save_pickle(rev_dict)

def save_pickle(rev_dict):
    with open(cache_pickle,'w') as f:
        pickle.dump(rev_dict,f)    

def mark_interesting(rev_dict,only_new=True):
    revisions=sorted(rev_dict.keys())
    for revision in revisions:
        rev=rev_dict[revision]
        if only_new and rev.interesting==None:
            print '\n'*50
            print 'REVISION:'
            print '%d %s'%(rev.rev_id,rev.author)

            # print link if in dict
            if rev.rev_id-1 in rev_dict:
                parent_hash=rev_dict[rev.rev_id-1].sha1
                print github_compare+parent_hash+'...'+rev.sha1
            if rev.status==None:
                print 'Test status: unavailable'
            elif rev.status=='':
                print 'Test status: Passed'
            else:
                print 'Test status: %s'%(rev.status)
            print rev.commit_message
        
            rev.interesting=input_yes_no('Is this an interesting commit?')
            rev_dict[revision]=rev
            save_pickle(rev_dict)

def update_cache(rev_dict,git_repository):
    r = Reporter('pulling revision data from test server RSS')

    feed = feedparser.parse(test_server_rss)

    new_revisions=0
    updated_revisions=0
    for item_count,item in enumerate(feed['items']):
        title = item['title'].split(':')
        author = title[0].strip()
        rev_id = long(title[1].strip())
        test_link = item['link']
        sha1 = get_hash_from_github_url( item['source']['href'] )

        # Get one letter status codes
        status=''
        summary_soup = BeautifulSoup(item['summary'])
        text =  summary_soup.get_text().split()
        status_index=0
        for i,x in enumerate(text):
            if x=='Status:':
                status_index=i
                break
        for l in text[status_index+1:]:
            # Break if not a one letter status code
            if len(l)!=1:
                break
            elif l in ['B','U','I','S']:
                status+=l
 
        # Create/update RosettaRevision object, if necessary
        if rev_id in rev_dict:
            revision_updated=False
            rev=rev_dict[rev_id]
            assert(rev.sha1==sha1)
            if rev.status!=status:
                rev.status=status
                revision_updated=True
            if revision_updated:
                updated_revisions+=1
                rev_dict[rev_id]=rev
        else:
            rev=RosettaRevision(rev_id)
            rev.status=status
            rev_dict[rev_id]=rev
            new_revisions+=1

        r.report(item_count)

    r.done()

    max_rev=first_git_rev_id
    for rev_id in rev_dict.keys():
        if rev_id>max_rev:
            max_rev=rev_id

    repo = Repo(git_repository)

    r = Reporter('checking for uncached revisions and updating information from repo')
    for rev_id in xrange(first_git_rev_id,max_rev+1):
        if rev_id not in rev_dict:
            rev=RosettaRevision(rev_id)
            rev_dict[rev_id]=rev
            new_revisions+=1

        rev=rev_dict[rev_id]
        # Check for attributes that were added in later versions of script
        if not hasattr(rev,'parents'):
            rev.parents=None
        if not hasattr(rev,'commit_message'):
            rev.commit_message=None

        # Update from git
        rev.get_info_from_repo(repo)

        rev_dict[rev_id]=rev
        r.report(rev_id)

    r.done()

    print 'Added to cache:\n  %d new revisions\n  %d updated revisions\n'%(new_revisions,updated_revisions)
    return rev_dict

def string_from_list(l):
    return_string=''
    for item in l:
        return_string+=item
    return return_string

def escape_for_latex(s):
    return_string=''
    for char in s:
        if char=='_':
            return_string+='\_'
        else:
            return_string+=char
    return return_string

def frame_from_revision(rev,parent_hash):
    # Helper function for output_slides
    time_string=rev.date_time.strftime('%Y-%m-%d %I:%M %p')

    frame_list=['\\begin{frame}[t,fragile,allowframebreaks]{Revision %d}\n\n'%(rev.rev_id)]

    frame_list.append('\\textbf{%s} --- {\\small %s --- \\href{%s}{%s}}\n\n'%(escape_for_latex(rev.author),time_string,github_compare+parent_hash+'...'+rev.sha1,rev.short_sha1))

    tests_passed=False
    if rev.status==None:
        test_string='Unavailable'
    elif rev.status=='':
        tests_passed=True
    else:
        test_string=''
        if 'B' in rev.status:
            test_string+='Build '
        if 'U' in rev.status:
            test_string+='Unit '
        if 'I' in rev.status:
            test_string+='Integration '
        if 'S' in rev.status:
            test_string+='Scorefunction '
        if test_string=='':
            test_string='Unknown'

    if tests_passed:
        frame_list.append('\\textit{\\href{%s%d}{All tests passed}}\n\n'%(test_server_url,rev.rev_id))
    else:
        frame_list.append('\\textit{\small \\href{%s%d}{Test changes: %s}}\n\n'%(test_server_url,rev.rev_id,test_string))

    lines=[s.strip() for s in rev.commit_message.split('\n')]
    main_lines=[]
    for i,line in enumerate(lines):
        if line.startswith('NEWS:'):
            news_line=line[5:].strip()
            if news_line=='':
                news_line=lines[i+1]
            frame_list.append('\\textbf{NEWS:} %s\n\n'%(escape_for_latex(news_line)))
        main_lines.append(line+'\n')

    frame_list.append('\\begin{lstlisting}\n%s\\end{lstlisting}\n'%(string_from_list(main_lines)))

    frame_list.append('\n\\end{frame}\n\n')
    return frame_list

def output_slides(file_location,rev_dict,revisions):
    if file_location.endswith('.pdf'):
        file_location=file_location[:-4]

    latex_list=['\\documentclass{beamer}\n\\usetheme{default}\n\\usepackage{hyperref}\n\\usepackage{listings}\n\\lstset{breaklines=true}\n\\lstset{basicstyle=\\tiny\\ttfamily}\n\\setbeamertemplate{frametitle continuation}{}\n\\setbeamertemplate{navigation symbols}{}%remove navigation symbols\n\\begin{document}\n']

    revisions.sort()
    r=Reporter('creating LaTeX for slides')
    for i,revision in enumerate(revisions):
        rev=rev_dict[revision]
        parent_hash=rev_dict[rev.rev_id-1].sha1
        latex_list.extend(frame_from_revision(rev,parent_hash))
        r.report(i)

    latex_list.append('\\end{document}')

    input_string=string_from_list(latex_list)

    # Run pdflatex
    latex_file=file_location+'.tex'
    with open(latex_file,'wb') as f:
        f.write(input_string.encode("UTF-8"))

    print 'Running pdflatex'
    p=subprocess.Popen(['pdflatex',latex_file,'-jobname',file_location],stdout=PIPE,stderr=PIPE)
    stdoutdata, stderrdata = p.communicate()

    print stderrdata

    r.done()

    print 'Saved pdf output'

def initialize_data():
    rev_dict={}
    return rev_dict

if __name__ == "__main__":
    main()
