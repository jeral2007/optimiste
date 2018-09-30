import subprocess as sp
from time import sleep
delay = 0.5
class Node(object):
    def __init__(self, node_dir, calc_files, strict=True):
        self.node_dir = node_dir
        self.calc_files = [node_dir+'/' + filename for filename in calc_files]
        for f, nf in zip(calc_files, self.calc_files):
            sp.call('mkdir {}'.format(node_dir), shell=True)
            sp.call('cp {} {}'.format(f, nf), shell=True)
        self.node_dir = node_dir
        self.task = None
        self.strict = strict
        self.occupied = None

    def save_files(self):
        for filename in self.calc_files:
            sp.call('cp {} {}0'.format(filename, filename), shell=True)

    def restore_files(self):
        for filename in self.calc_files:
            sp.call('cp {}0 {}'.format(filename, filename), shell=True)

    def free(self):
        sp.call('rm -rf {}'.format(self.node_dir), shell=True)

    def run(self, command):
        if self.task is None or self.task.poll() == 0:
            self.task = sp.Popen(command, shell=True, cwd=self.node_dir)
        sleep(delay)

    def status(self):
        if self.task is None:
            return 'IDLE'
        elif self.task.poll() is None:
            return 'RUNNING'
        elif self.task.poll() == 0:
            return 'FINISHED'
        elif self.strict:
            return 'ERROR'
        else:
            return 'IDLE'


class Task(object):
    last_id = 0
    def __init__(self, node, commands):
        self.node = node
        self.result = [[]]
        self.commands = commands
        self.state = 0
        self.finished = False
        self.id = type(self).last_id
        type(self).last_id += 1

    def run(self):
        if self.finished:
            return
        if self.node.occupied is not None and self.node.occupied != self.id:
            return

        print "task {}, state = {}\
            last non-terminal state = {}".format(self.id,
                                                 self.state,
                                                 len(self.commands)-1)
        if  self.state >= len(self.commands):
            self.finished = True
            self.node.occupied = None
            return

        if not self.node.status() in ['IDLE', 'FINISHED']:
            return

        self.node.occupied = self.id
        command, blocking = self.commands[self.state]
        while blocking and self.state<len(self.commands)-1:
            command(self.node, self.result)
            self.state += 1
            command, blocking = self.commands[self.state]
        command(self.node, self.result)
        self.state +=1
        return
