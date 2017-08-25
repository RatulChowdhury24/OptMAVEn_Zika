""" This module reports the times that each step takes, to assess the
performance. """

import datetime
import os
import re
import sys
import tempfile
import time

import EXPERIMENT


class File(object):
    """ A file that records the performance of an OptMAVEn 2.0
    experiment. """
    def __init__(self, filename, new_file=False):
        self.filename = filename
        # Make a new file if instructed or if the file does not exist.
        if new_file or not os.path.isfile(filename):
            open(filename, "w")
        self.start_key = "START"
        self.finish_key = "FINISH"
        self.time_key = "TIME"
        self.sub_time_key = "SUB_TIME"
        self.du_key = "DISK_USAGE"
        self.aspects = (self.start_key, self.finish_key, self.time_key,
                self.sub_time_key, self.du_key)
        self.dtime_format = "on %Y %B %d at %H:%M:%S" 

    def read(self):
        content = open(self.filename).read()
        return content

    def write(self, aspect, process, value):
        # Record an aspect of the performance of the experiment.
        if aspect not in self.aspects:
            raise KeyError("Bad aspect for performance file: {}".format(aspect))
        open(self.filename, "a").write("{}\t{}\t{}\n".format(aspect, process,
                value))
    
    def start(self, experiment, t=None):
        if self.read().strip() != "":
            raise ValueError("Cannot start experiment that is already started.")
        if not isinstance(t, datetime.datetime):
            t = datetime.datetime.now()
        self.write(self.start_key, experiment["Name"], t.strftime(
                self.dtime_format))
    
    def finish(self, experiment, t=None):
        if not isinstance(t, datetime.datetime):
            t = datetime.datetime.now()
        self.write(self.finish_key, experiment["Name"], t.strftime(
                self.dtime_format))

    def time(self, process, time_):
        if isinstance(time_, Timer):
            time_ = time_.time()
        else:
            # time_ should be in seconds.
            time_ = float(time_)
        self.write(self.time_key, process, time_)

    def sub_time(self, process, time_):
        if isinstance(time_, Timer):
            time_ = time_.time()
        else:
            # time_ should be in seconds.
            time_ = float(time_)
        self.write(self.sub_time_key, process, time_)

    def disk_usage(self, process, directory):
        handle, du_file = tempfile.mkstemp()
        os.system("du -sk {} > {}".format(directory, du_file))
        du = int(open(du_file).read().split()[0])
        try:
            os.remove(du_file)
        except OSError:
            pass
        self.write(self.du_key, process, du)        
    
    def duration(self):
        """ Time difference between finish and start times. """
        start = None
        finish = None
        split = self.dtime_format.split()[0]
        with open(self.filename) as f:
            for line in f:
                if line.startswith((self.start_key, self.finish_key)):
                    stamp = (split + line.split(split)[-1]).strip("\n")
                    dtime = datetime.datetime.strptime(stamp, self.dtime_format)
                    if line.startswith(self.start_key):
                        if start is not None:
                            raise ValueError("Cannot start twice")
                        start = dtime
                    else:
                        if finish is not None:
                            raise ValueError("Cannot finish twice")
                        finish = dtime
        return (finish - start).total_seconds()
    
    def clock_time(self):
        """ Sum of all CPU times (user + sys, using Unix time utility). """
        return sum(map(float, [line.split("\t")[2] for line in open(
                self.filename) if line.startswith(self.time_key)]))

    def get_sub_time(self, process, allow_duplicates=False,
            error_if_not_found=True):
        sub_time = None
        with open(self.filename) as f:
            for line in f:
                if line.startswith(self.sub_time_key):
                    if process in line:
                        t = float(line.split(process)[1])
                        if sub_time is None:
                            sub_time = t
                        elif allow_duplicates:
                            sub_time += t
                        else:
                            raise ValueError("Duplicate of {}".format(process))
        if sub_time is None:
            if error_if_not_found:
                raise ValueError("Cannot find {}".format(process))
            else:
                sub_time = 0.0
        return sub_time

    def time_process_category(self, category, count="all", indexes=None):
        if count not in ["all", "sub_time", "time"]:
            raise ValueError("Bad count argument: {}".format(count))
        time_ = 0.0
        starts = tuple([self.sub_time_key] * int(count in ["all", "sub_time"]) +
                [self.time_key] * int(count in ["all", "time"]))
        with open(self.filename) as f:
            index = 0
            for line in f:
                if line.startswith(starts):
                    if category in line:
                        if indexes is None or index in indexes:
                            time_ += float(line.split()[-1])
                        index += 1
        return time_
 

class Timer(object):
    def __init__(self, report_file=None):
        self.begin = None
        self.end = None
        self.running = False
        self.what = None
        self.file = report_file
    
    def start(self, what=None):
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
    
    """ 
    def dump(self, report_file=None, mode="a"):
        if self.file is None:
            rprt = report_file
        else:
            rprt = self.file
        if rprt is None:
            raise ValueError("Must have a file to report time.")
        open(rprt, mode).write("\n{} took {} seconds.\n".format(self.what,
                self.time()))
    """

def performance_file(experiment):
    """ Return the performance file of the experiment. """
    return File(os.path.join(experiment["Folder"], "performance.txt"))


if __name__ == "__main__":
    time_file = sys.argv[1]
    process = sys.argv[2]
    time_ = 0.0
    with open(time_file) as f:
        for line in f:
            if line.startswith(("user", "sys")):
                if len(line.split()) == 2:
                    try:
                        time_ += float(line.split()[1])
                    except ValueError:
                        pass
    p_file = performance_file(EXPERIMENT.Experiment())
    p_file.time(process, time_)
