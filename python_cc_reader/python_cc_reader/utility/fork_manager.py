# basically, I'm borrowing/stealling some of Sergey's awesome logic
# for controlling forked processes, and putting it into the
# python_cc_reader directory
import os, sys
import time

class ForkManager:
    def __init__( self, max_jobs = 1 ) :
        self.jobs = []
        self.max_n_jobs = max_jobs
        self.success_callback = None
        self.error_callback = None

    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
        '''
        while len(self.jobs) >= self.max_n_jobs :
            for p in self.jobs[:] :
                r = os.waitpid(p, os.WNOHANG)
                if r == (p, 0):  # process has ended without error
                    self.jobs.remove(p)
                    if self.success_callback : self.success_callback( self, p )
                elif r[0] == p :  # process ended but with an error
                    self.jobs.remove(p)
                    if self.error_callback : self.error_callback( self, p )

            if len(self.jobs) >= self.max_n_jobs : time.sleep(.1)

        sys.stdout.flush()
        pid = os.fork()
        if pid: self.jobs.append(pid) # We are parent!
        return pid

    def wait_for_remaining_jobs( self ) :
        ''' Block until all child processes have finished
        '''
        for p in self.jobs:
            r = os.waitpid(p, 0)
            if r == (p,0) : # process has ended without error
                if self.success_callback : self.success_callback( self, p )
            elif r[0] == p : # process has ended with an error
                if self.error_callback : self.error_callback( self, p )
