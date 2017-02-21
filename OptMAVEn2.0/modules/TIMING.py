""" This module reports the times that each step takes, to assess the
performance. """

import time


class Timer(object):
    def __init__(self, report_file=None):
        self.begin = None
        self.end = None
        self.running = False
        self.what = None
        self.file = report_file
    
    def start(self, what):
        self.end = None
        self.what = what
        self.begin = time.time()
        self.running = True
    
    def stop(self):
        if not self.running:
            raise ValueError("Cannot stop timer that is not running.")
        self.end = time.time()
        self.running = False
        return self.time()
    
    def time(self):
        if self.begin is not None:
            if self.running:
                return time.time() - self.begin
            else:
            	return self.end - self.begin
        else:
            raise ValueError("The timer was not started.")
    
    def dump(self, report_file=None, mode="a"):
        if self.file is None:
            rprt = report_file
        else:
            rprt = self.file
        if rprt is None:
            raise ValueError("Must have a file to report time.")
        open(rprt, mode).write("{} took {} seconds.\n".format(self.what,
                self.time()))
