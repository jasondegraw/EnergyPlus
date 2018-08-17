import os
import re
import filevisitor

def mergetxt(oldtext, oldlabel, newtext, newlabel):
    return '<<<<<<< %s\n%s\n=======\n%s\n>>>>>>> %s' % (oldlabel,
                                                        oldtext,
                                                        newtext,
                                                        newlabel)

class ArgGetter:
    def __init__(self):
        self.complete = False
        self.inquote = False
        self.inchar = False
        self.string = ''
    def start(self, line):
        line = line.partition('(')[2]
        self.count = 1
        self.resume(line)
    def resume(self, line):
        it = iter(line)
        for c in it:
            if self.inquote:
                if c == '\\': # Skip escaped stuff
                    self.string += c
                    c = next(it)
                elif c == '"':
                    self.inquote = False
            elif self.inchar:
                if c == '\\': # Skip escaped stuff
                    self.string += c
                    c = next(it)
                elif c == "'":
                    self.inchar = False
            else:
                if c == '"':
                    self.inquote = True
                elif c == "'":
                    self.inchar = True
                elif c == ')':
                    self.count -= 1
                    if self.count == 0:
                        self.complete = True
                        return;
                elif c == '(':
                    self.count += 1
            self.string += c
        self.string += '\n'
    def secondresume(self, line):
        it = iter(line)
        for c in it:
            if c == '\\': # Skip escaped stuff
                self.string += c
                c = next(it)
            elif c == '"' or c == "'":
                self.inquote = not self.inquote
            elif c == ')':
                if not self.inquote:
                    self.count -= 1
                    if self.count == 0:
                        self.complete = True
                        return;
            elif c == '(':
                if not self.inquote:
                    self.count += 1
            self.string += c
        self.string += '\n'
    def firstresume(self, line):
        for c in line:
            if c == ')':
                self.count -= 1
                if self.count == 0:
                    self.complete = True
                    return;
            elif c == '(':
                self.count += 1
            self.string += c
        self.string += '\n'

class FileChecker(filevisitor.CodeChecker):
    def __init__(self, string):
        filevisitor.CodeChecker.__init__(self)
        self.string = string
        self.files = []
        self.filehits = {}
    #def filter(self, filepath):
    #    if 'UtilityRoutines.' in filepath:
    #        return False
    #    return filepath.endswith('.cc') or filepath.endswith('.hh')
    def visit(self, filepath):
        fp = open(filepath, 'r')
        hits = 0
        for line in fp:
            if line.lstrip().startswith('//'):
                continue
            if self.string in line:
                hits += 1
        #text = fp.read()
        #fp.close()
        #hits = text.count(self.string)
        if hits > 0:
            self.files.append(filepath)
            self.filehits[filepath] = hits
        return hits

class EmptyChecker(FileChecker):
    def __init__(self, string):
        FileChecker.__init__(self, string)
        self.empty = 0
    def visit(self, filepath):
        fp = open(filepath, 'r')
        hits = 0
        for line in fp:
            if line.lstrip().startswith('//'):
                continue
            if self.string in line:
                hits += 1
                # Grab the arguments
                getter = ArgGetter()
                getter.start(line)
                while not getter.complete:
                    getter.resume(next(fp))
                #print(filepath, getter.string)
                print(getter.string)
                if getter.string.strip() == '""':
                    self.empty += 1
        if hits > 0:
            self.files.append(filepath)
        return hits

class UnexpectedCondition(Exception):
    pass

class ErrorCall:
    def __init__(self, message, line):
        self.message = message
        self.line = line
    @classmethod
    def type(cls):
        return cls.__name__
    def __str__(self):
        return self.type()+'('+self.message+');'

class StartError(ErrorCall):
    def __init__(self, message, line):
        ErrorCall.__init__(self, message, line)
        self.continuations = []
    def append(self, continuation):
        self.continuations.append(continuation)
    def __str__(self):
        text = ErrorCall.__str__(self)
        for el in self.continuations:
            text += '\n' + str(el)
        return text
    def pattern(self):
        lines = str(self).splitlines()
        #lines.extend([str(el) for el in self.continuations])
        return '\s*' + '\s*'.join([re.escape(el) for el in lines])
        
class ContinueError(ErrorCall):
    pass

class ShowFatalError(StartError):
    def __init__(self, message, line):
        super().__init__(message, line)
    def reformat(self, indent=''):
        text = 'fatal( %s );' % self.message.strip()
        lines = [el.lstrip() for el in text.splitlines()]
        hanging = '\n' + indent + '\t'
        return indent + hanging.join(lines)

class ShowSevereError(StartError):
    def __init__(self, message, line):
        super().__init__(message, line)
    def reformat(self, indent=''):
        fcn = 'error'
        if self.continuations:
            vec = []
            timestampMsg = ''
            for cont in self.continuations:
                if cont.type() == 'ShowContinueErrorTimeStamp':
                    if cont.message.strip() != '""':
                        if timestampMsg:
                            raise UnexpectedCondition("Duplicate time stamp calls")
                        timestampMsg = cont.message.strip()
                else:
                    vec.append(cont.message.strip())
            #vec = [el.message.strip() for el in self.continuations]
            #print(vec)
            for sub in self.continuations:
                if sub.type() == 'ShowFatalError': # Upgrade to fatal
                    fcn = 'fatal'
                    break
            if timestampMsg:
                timestampMsg = ',\n' + timestampMsg
            text = '%s( %s,\n{ %s }%s );' % (fcn,
                                             self.message.strip(),
                                             ',\n'.join(vec),
                                             timestampMsg)
        else:
            text = '%s( %s );' % (fcn, self.message.strip())
        lines = [el.strip() for el in text.splitlines()]
        #print(lines)
        hanging = '\n' + indent + '\t'
        return indent + hanging.join(lines)

class ShowSevereMessage(StartError):
    pass

class ShowContinueError(ContinueError):
    pass

class ShowContinueErrorTimeStamp(ContinueError):
    pass

class ShowMessage(StartError):
    pass

class ShowWarningError(StartError):
    def __init__(self, message, line):
        super().__init__(message, line)
    def reformat(self, indent=''):
        fcn = 'warning'
        if self.continuations:
            vec = []
            timestampMsg = ''
            for cont in self.continuations:
                if cont.type() == 'ShowContinueErrorTimeStamp':
                    #print('####%s####' % cont.message)
                    if cont.message.strip() != '""':
                        if timestampMsg:
                            raise UnexpectedCondition("Duplicate time stamp calls")
                        timestampMsg = cont.message.strip()
                else:
                    vec.append(cont.message.strip())
            #vec = [el.message.strip() for el in self.continuations]
            if timestampMsg:
                timestampMsg = ',\n' + timestampMsg
            if vec:
                text = '%s( %s,\n{ %s }%s );' % (fcn,
                                                 self.message.strip(),
                                                 ',\n'.join(vec),
                                                 timestampMsg)
            else:
                text = '%s( %s%s );' % (fcn,self.message.strip(),timestampMsg)
        else:
            text = '%s( %s );' % (fcn, self.message.strip())
        lines = [el.lstrip() for el in text.splitlines()]
        hanging = '\n' + indent + '\t'
        return indent + hanging.join(lines)

class ShowWarningMessage(StartError):
    pass

class ShowRecurringSevereErrorAtEnd(StartError):
    pass

class ShowRecurringWarningErrorAtEnd(StartError):
    pass

class ShowRecurringContinueErrorAtEnd(ContinueError):
    pass

class ShowErrorMessage(StartError):
    pass

class ShowRecurringErrors(StartError):
    pass    

class ErrorFinderName(filevisitor.CodeChecker):
    def __init__(self, starter, continuations):
        filevisitor.CodeChecker.__init__(self)
        self.starter = starter
        self.files = []
        self.continuations = continuations
        self.continued = 0
        self.indeterminate = 0
        self.calls = {}
    #def filter(self, filename):
    #    return filename.endswith('AirflowNetworkBalanceManager.cc')
    #    #return filename.endswith('ConvectionCoefficients.cc')
    def visit(self, filepath):
        fp = open(filepath, 'r')
        lineNumber = 0
        calls = []
        hits = 0
        for line in fp:
            lineNumber += 1
            if line.lstrip().startswith('//'):
                continue
            if self.starter in line:
                hits += 1
                # Grab the arguments
                print(self.starter, lineNumber)
                getter = ArgGetter()
                getter.start(line)
                while not getter.complete:
                    lineNumber += 1
                    getter.resume(next(fp))
                calls.append(self.starter+'('+getter.string+');')
                # Look for continuations
                lineNumber += 1
                nextLine = next(fp)
                if self.starter in nextLine:
                    raise UnexpectedCondition('Unexpected error starter call')
                multiLine = False
                notDone = True
                while notDone: # Yay!
                    for cont in self.continuations:
                        if nextLine.lstrip().startswith(cont + '('):
                            print(cont)
                            cogetter = ArgGetter()
                            cogetter.start(nextLine)
                            while not cogetter.complete:
                                cogetter.resume(next(fp))
                            calls[-1] += '\n' + cont+'('+cogetter.string+');'
                            multiLine = True
                            lineNumber += 1
                            nextLine = next(fp)
                            break
                    else:
                        notDone = False
                if multiLine:
                    self.continued += 1
        self.calls[filepath] = calls
        #text = fp.read()
        #fp.close()
        #hits = text.count(self.string)
        if hits > 0:
            self.files.append(filepath)
        return hits

class ErrorFinder(filevisitor.CodeChecker):
    def __init__(self, starter, continuations):
        filevisitor.CodeChecker.__init__(self)
        self.starter = starter
        self.files = []
        self.filehits = {}
        self.continuations = continuations
        self.continued = 0
        self.indeterminate = 0
        self.calls = {}
    #def filter(self, filename):
    #    return filename.endswith('AirflowNetworkBalanceManager.cc')
    #    #return filename.endswith('ConvectionCoefficients.cc')
    def visit(self, filepath):
        fp = open(filepath, 'r')
        lineNumber = 0
        calls = []
        hits = 0
        for line in fp:
            lineNumber += 1
            if line.lstrip().startswith('//'):
                continue
            if self.starter.type() in line:
                hits += 1
                startNumber = lineNumber
                # Grab the arguments
                getter = ArgGetter()
                getter.start(line)
                while not getter.complete:
                    getter.resume(next(fp))
                    lineNumber += 1
                calls.append(self.starter(getter.string,
                                          startNumber))
                # Look for continuations
                nextLine = next(fp)
                lineNumber += 1
                #if self.starter.type() in nextLine:
                #    raise UnexpectedCondition('Unexpected error starter call')
                multiLine = False
                notDone = True
                while notDone: # Yay!
                    for cont in self.continuations:
                        if nextLine.lstrip().startswith(cont.type() + '('):
                            contNumber = lineNumber
                            cogetter = ArgGetter()
                            cogetter.start(nextLine)
                            while not cogetter.complete:
                                cogetter.resume(next(fp))
                                lineNumber += 1
                            calls[-1].append(cont(cogetter.string,
                                                  contNumber))
                            #if cont.type() == 'ShowContinueErrorTimeStamp':
                            #    print('$$$$$$$')
                            multiLine = True
                            nextLine = next(fp)
                            lineNumber += 1
                            break
                    else:
                        notDone = False
                if multiLine:
                    self.continued += 1
        self.calls[filepath] = calls; 
        #text = fp.read()
        #fp.close()
        #hits = text.count(self.string)
        if hits > 0:
            self.files.append(filepath)
            self.filehits[filepath] = hits
        return hits

allnames = '''ShowFatalError
ShowSevereError
ShowSevereMessage
ShowContinueError
ShowContinueErrorTimeStamp
ShowMessage
ShowWarningError
ShowWarningMessage
ShowRecurringSevereErrorAtEnd
ShowRecurringWarningErrorAtEnd
ShowRecurringContinueErrorAtEnd
StoreRecurringErrorMessage
ShowErrorMessage
SummarizeErrors
ShowRecurringErrors'''.splitlines()

objnames = '''ShowFatalError
ShowSevereError
ShowSevereMessage
ShowMessage
ShowWarningError
ShowWarningMessage
ShowErrorMessage
ShowContinueError
ShowContinueErrorTimeStamp'''.splitlines()

path = '../../src/EnergyPlus/'
filepath = '../../src/EnergyPlus/AirflowNetworkBalanceManager.cc'

lookup = {}
for fcn in objnames:
    lookup[fcn] = globals()[fcn]

fp = open('counts.csv', 'w')

# Do a basic search to get usage counts
fp.write('Name, Usage, Source Files, Header Files\n')
for fcn in allnames:
    checker = FileChecker(fcn)
    checker.execute(path)
    hh = 0
    for file in checker.filehits.keys():
        if file.endswith('.hh'):
            hh += 1
    fp.write('%s,%d,%d,%d\n' % (fcn, checker.count, len(checker.files)-hh, hh))
    #if fcn == 'ShowWarningMessage':
    #    for file, count in checker.filehits.items():
    #        print(file, count)
fp.close()

startObj = [ShowSevereError,
            ShowSevereMessage,
            ShowMessage,
            ShowWarningError,
            ShowWarningMessage,
            ShowFatalError]

continueObj = [ShowContinueError,
               ShowContinueErrorTimeStamp,
               ShowFatalError]

finders = []

for starter in startObj:
    print('Looking for %s continued errors...' % starter.type())
    finder = ErrorFinder(starter, continueObj)
    #finder.execute(path)
    finder.execute_once(filepath)
    finders.append(finder)
    print('\tFound %d' % finder.count)

##fp = open('results.csv','w')
##
##for fcn in allnames.splitlines():
##    checker = FileChecker(fcn)
##    checker.execute(path)
##    hh = 0
##    for file in checker.files:
##        if file.endswith('.hh'):
##            hh += 1
##    fp.write('%s,%d,%d,%d\n' % (fcn, checker.count, len(checker.files)-hh, hh))
##    print('%s,%d,%d,%d' % (fcn, checker.count, len(checker.files)-hh, hh))
##    if fcn in headercalls:
##        for file in checker.files:
##            if file.endswith('.hh'):
##                print('\t',file)
##    if fcn in lowcalls:
##        for file in checker.files:
##            print('\t',file)
##
##fp.close()

##finder = ErrorFinderName('ShowSevereError', ['ShowContinueError',
##                                             'ShowContinueErrorTimeStamp'])
##finder.execute(path)
##print(finder.count, finder.continued)
##
##findErrors = ErrorFinder(ShowSevereError, [ShowContinueError,
##                                           ShowContinueErrorTimeStamp,
##                                           ShowFatalError])
##findErrors.execute(path)
##print(findErrors.count, findErrors.continued)

##findWarnings = ErrorFinder(ShowWarningError, [ShowContinueError,
##                                              ShowContinueErrorTimeStamp,
##                                              ShowFatalError])
##findWarnings.execute(path)
##print(findWarnings.count, findWarnings.continued)
##
### Combine everything into one big bucket
##keys = set(list(findErrors.calls.keys()) + list(findWarnings.calls.keys()))
##
##calls = {}
##for key in keys:
##    calls[key] = findErrors.calls[key] + findWarnings.calls[key]
##
##for filepath, calls in finder.calls.items():
##    print(filepath, len(calls))
##    fp = open(filepath, 'r')
##    text = fp.read()
##    fp.close()
##    count = 0
##    matched = 0
##    for call in calls:
##        lines = call.splitlines()
##        if len(lines) > 1:
##            rex = '\s*' + '\s*'.join([re.escape(el) for el in lines])
##            #print(rex)
##            #print()
##            m = re.search(rex,text)
##            if m:
##                #if 'TimeStamp' in m.group(0):
##                #    print(m.group(0))
##                #print(m.group(0))
##                matchtxt = m.group(0)
##                occ = text.count(matchtxt)
##                #if occ > 1:
##                #    print(matchtxt)
##                #    raise UnexpectedCondition('Multiple instances (%d)!' % occ)
##                matched += 1
##            else:
##                print(call)
##        elif call in text:
##            #print(call)
##            count += 1
##        else:
##            print(call)
##    print(filepath, count, len(calls))
##
##for filepath, calls in findErrors.calls.items():
##    print(filepath, len(calls))
##    fp = open(filepath, 'r')
##    text = fp.read()
##    fp.close()
##    count = 0
##    matched = 0
##    duplicate = 0
##    for start in calls:
##        rex = start.pattern()
##        m = re.search(rex,text)
##        if m:
##            matchtxt = m.group(0)[1:]
##            occ = text.count(matchtxt)
##            if occ > 1:
##                duplicate += 1
##                #print(matchtxt)
##                #print('Multiple instances (%d)!' % occ)
##                #raise UnexpectedCondition('Multiple instances (%d)!' % occ)
##            else:
##                indent = ''
##                for c in matchtxt:
##                    if c.isspace():
##                        indent += c
##                    else:
##                        break
##                rep = mergetxt(matchtxt, 'old', start.reformat(indent), 'new')
##                text = text.replace(matchtxt, rep)
##            matched += 1
##        else:
##            print('++++')
##            print(start)
##    print(filepath, matched, len(calls))
##    #fp = open(filepath, 'w')
##    #fp.write(text)
##    #fp.close()
##
###empty = EmptyChecker('ShowContinueErrorTimeStamp')
###empty.execute(path)
##
##filepath = path + 'AirflowNetworkBalanceManager.cc'
##
##for i, call in enumerate(findErrors.calls[filepath]):
##    print(i, ',', call.type(), ',', call.line)
##
###print(findex.calls[filepath][0].reformat('\t'))
