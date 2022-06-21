#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Based on: https://selvakumar-arumugapandian.medium.com/command-line-subcommands-with-pythons-argparse-4dbac80f7110

import os
import sys
import shutil
import argparse
from pathlib import Path

import pynteny.subcommands as sub



class Pynteny():

    def __init__(self):
        parser = argparse.ArgumentParser(
            description=("Pynteny: synteny-based hmm searches made easy in Python"),
            epilog="Semidán Robaina Estévez (srobaina@ull.edu.es), 2022"
            )
        parser.add_argument("command", help="pynteny subcommand")
        parser.add_argument("-v","--version", help="show version and exit", action="version", version="0.0.1")
        
        #Read the first argument
        args = parser.parse_args(sys.argv[1:2])
        #use dispatch pattern to invoke method with same name of the argument
        getattr(self, args.command)()
    
    def seach(self):
        parser = argparse.ArgumentParser(description="Adds a file")
        parser.add_argument()

        #we are inside a subcommand, so ignore the first argument and read the rest
        args = parser.parse_args(sys.argv[2:])
        sub.synteny_search(args)
    
    def preprocess(self):
        parser = argparse.ArgumentParser(description="Commits a file")
        parser.add_argument()
        #we are inside a subcommand, so ignore the first argument and read the rest
        args = parser.parse_args(sys.argv[2:])
        sub.translate_assembly(args)


