#!/usr/bin/env python
###############################################################################
### imports
###############################################################################
import sys
import os
import subprocess
import time
import datetime

###############################################################################
### utility functions
###############################################################################
def usage():
    print "usage:", __file__

def format_time(t):
    t = datetime.datetime.fromtimestamp(t)
    return '[{time}]'.format(time=t.strftime('%Y-%m-%d %H:%M:%S'))

def log(*m, **kw):
    m = list(m)
    if 'time' in kw:
        m.insert(0, format_time(kw['time']))
    sys.stdout.write(' '.join(m))
    sys.stdout.write('\n')
    sys.stdout.flush()
    return

def check_output(*args, **kwargs):
    for k, v in kwargs.iteritems():
        args.append(' '.join(k, v))
    o, e = subprocess.Popen(
        ' '.join(args), shell=True, executable='/bin/bash',
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    ).communicate()
    if e and len(e):
        #raise Exception(e)
        print e
    return o


###############################################################################
### helper classes
###############################################################################
class AutoQueue(object):

    def __init__(self, 
                 status_cmd = None, 
                 submit_cmd = None, 
                 submit_file = None
    ):
        self._user = os.getlogin()
        self._status_cmd = status_cmd.replace('USER', self._user)
        self._submit_cmd = submit_cmd.replace('FILE', submit_file)
        self._submit_file = submit_file
        if not os.path.exists(submit_file):
            raise Exception(submit_file+" does not exist!!!")
        self._sleep_time = 1*60
        self._queue = []

    def _submit_jobs(self):
        '''
        details: source submit_file
        
        '''
        o = check_output(self._submit_cmd)
        self._queue = filter(None, o.split('\n'))

    def _update_queue(self):
        '''
        details: run status command and parse jobs in queue
        
        '''
        o = check_output(self._status_cmd)
        self._queue = [id for id in self._queue if id in o]

    def perpetuate(self):
        '''
        details: this method watches individual job ids, and re-submits these 
                 jobs whenever they are no longer found in the queue
 
        '''
        while True:
            try:
                self._update_queue()
                if len(self._queue):
                    time.sleep(self._sleep_time)
                    continue
                log('queue is empty', time=time.time())
                self._submit_jobs()
                log(self._submit_file,
                    'has been submitted\n',
                    time=time.time())
            except Exception as e:
                log(e.__str__())
                time.sleep(self._sleep_time)
        return False

    def _cluster_queue_is_empty(self):
        '''
        details: run status_cmd and check output
        
        '''
        o = check_output(self._status_cmd)
        o = filter(None, o.split('\n'))
        return (len(o) < 1)

    def perpetuate_nonspecific(self):
        '''
        details: this method watches the user's entire queue, and re-submits  
                 jobs whenever it is empty (i.e. no jobs returned by qstat)
 
        todo: test, deprecate?
        '''
        while True:
            try:
                if not self._cluster_queue_is_empty():
                    time.sleep(self._sleep_time)
                    continue
                log('cluster queue is empty', time=time.time())
                self._submit_jobs()
                log(self._submit_file,
                    'has been submitted\n',
                    time=time.time())
            except Exception as e:
                log(e.__str__())
                time.sleep(self._sleep_time)
        return False


class SLURMAutoQueue(AutoQueue):

    def __init__(self):
        AutoQueue.__init__(
            self,
            status_cmd = "squeue --user USER | tail --lines=+2",
            submit_cmd = "source FILE | awk '{print $4}'",
            submit_file = "./sbatchMINI"
        )


class PBSAutoQueue(AutoQueue):

    def __init__(self):
        AutoQueue.__init__(
            self,
            status_cmd = "qstat | tail --lines=+3",
            submit_cmd = "source FILE | awk '{print $1}'",
            submit_file = "./qsubMINI"
        )


###############################################################################
### factories
###############################################################################
class AutoQueueFactoryException(Exception):
    pass
    
class AutoQueueFactory(object):
    
    def __init__(self):
        pass
    
    def get_auto_queue(self):
        if os.path.exists('/usr/bin/sbatch'):
            return SLURMAutoQueue()
        elif os.path.exists('/usr/bin/qsub'):
            return PBSAutoQueue()
        else:
            raise AutoQueueFactoryException("Could not find sbatch or qsub!!!")

auto_queue_factory = AutoQueueFactory()


###############################################################################
### main function
###############################################################################
def auto_queue(*args, **kwargs):

    ### TODO: should this just be a method of Queue?
    q = auto_queue_factory.get_auto_queue()
    q.perpetuate()


###############################################################################
### main script
###############################################################################
if __name__ == '__main__':

    args = sys.argv[1:]
    kwargs = {}

    sys.exit(auto_queue(*args, **kwargs))
