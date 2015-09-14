#!/bin/bash

mv *.log $1/
mv *.root $1/
rm -r input
rm -r GRP
rm gphy*
rm QAImageRec*
rm validation_error.message
