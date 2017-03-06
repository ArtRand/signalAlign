#!/usr/bin/env python2.7

import os


class FolderHandler(object):
    def __init__(self):
        self.path = None  # place to make temp folder
        self.open = False  # flag for open or closed
        self.filenames = []  # place to hold paths to temp files

    def open_folder(self, path):
        self.path = os.path.expanduser(path)

        if not self.path.endswith("/"):
            self.path = self.path + "/"

        # if the folder already exists, make sure wr keep track so that we don't remove anything we
        # didn't want to
        if not os.path.isdir(path):
            os.system("mkdir {dir}".format(dir=self.path))
            self.already_exists = False
        else:
            self.already_exists = True

        self.open = True
        return self.path

    def add_file_path(self, filename):
        if self.open is True:
            if (self.path + filename) not in self.filenames:
                self.filenames.append(self.path + filename)
            return self.path + filename

    def remove_file(self, filename):
        if self.open is True:
            if (filename in self.filenames) and (os.path.isfile(filename)):
                os.remove(filename)
                self.filenames.remove(filename)
                return True
        else:
            return False

    def remove_folder(self):
        for _ in self.filenames:
            try:
                os.remove(_)
            except OSError:
                continue
        if self.already_exists:
            self.open = False
            return
        if os.listdir(self.path) == []:
            os.removedirs(self.path)
            self.open = False
            return


