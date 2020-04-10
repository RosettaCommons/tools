# basically, I'm borrowing/stealling some of Sergey's awesome logic
# for controlling forked processes, and putting it into the
# python_cc_reader directory
import os, sys
import time

class ForkManager:
    def __init__( self, max_jobs = 1 ) :
        self.jobs = [None] * max_jobs
        self.max_n_jobs = max_jobs
        self.success_callback = None
        self.error_callback = None
        self.myjob_index = None
        
    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
        '''
        avail_jobid = None
        while avail_jobid is None:
            for i, process_id in enumerate(self.jobs):
                if process_id is None:
                    avail_jobid = i
                    break
                try:
                    r = os.waitpid(process_id, os.WNOHANG)
                    if r[0] == process_id:
                        avail_job = i
                        self.jobs[i] = None
                        if r[1] == 0:
                            if self.success_callback:
                                self.success_callback(self, process_id)
                        else:
                            if self.error_callback:
                                self.error_callback(self, process_id)
                        break
                except OSError:
                    print("Warning: child process %d no found" % p)
                    avail_job = i
                    self.jobs[i] = None
                    if self.error_callback:
                        self.error_callback(self, process_id)
                    break
                
        # while len(self.jobs) >= self.max_n_jobs :
        #     for p in self.jobs[:] :
        #         try :
        #             r = os.waitpid(p, os.WNOHANG)
        #             if r == (p, 0):  # process has ended without error
        #                 self.jobs.remove(p)
        #                 if self.success_callback : self.success_callback( self, p )
        #             elif r[0] == p :  # process ended but with an error
        #                 self.jobs.remove(p)
        #                 if self.error_callback : self.error_callback( self, p )
        #         except OSError:
        #             print("Warning: child process %d not found" % p)
        #             self.jobs.remove(p)
        #             if self.error_callback : self.error_callback(self, p)
        # 
        #     if len(self.jobs) >= self.max_n_jobs : time.sleep(.1)

        # flush stdout before forking or output that has been buffered
        # but not delivered will be written by both the parent and
        # child processes!
        sys.stdout.flush()

        pid = os.fork()
        if pid:
            # We are the parent process!
            self.jobs[avail_jobid] = pid
        else:
            # we are the child process!
            self.myjob_index = avail_jobid
        return pid

    def wait_for_remaining_jobs( self ) :
        ''' Block until all child processes have finished
        '''
        for i, p in enumerate(self.jobs):
            if p is None:
                continue
            try:
                r = os.waitpid(p, 0)
                self.jobs[i] = None
                assert r[0] == p
                if r[1] == 0 : # process has ended without error
                    if self.success_callback :
                        self.success_callback( self, p )
                else:
                    # process has ended with an error
                    if self.error_callback : self.error_callback( self, p )
            except OSError:
                print("Warning: child process %d not found" % p)
                self.jobs[i] = None
                if self.error_callback:
                    self.error_callback(self, p)
