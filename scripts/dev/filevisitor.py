import os

class FileVisitor:
    def __init__(self):
        self.count = 0
    def visit(self, filepath):
        return 1
    def filter(self, filename):
        return true
    def files(self, path):
        return glob.glob(path+'*')
    def execute_once(self, file_path):
        if self.filter(file_path):
            self.count += self.visit(file_path)
    def execute(self, path, recursive=False):
        for path, subdirs, files in os.walk(path):
            for file in files:
                fullname = os.path.join(path,file)
                if self.filter(fullname):
                    self.count += self.visit(fullname)
            if not recursive:
                break

class CodeChecker(FileVisitor):
    def __init__(self):
        FileVisitor.__init__(self)
    def filter(self, filename):
        return filename.endswith('.cc') or filename.endswith('.hh')
