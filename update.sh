#!/bin/bash
cp -R ../MAPGD-DOXYGEN/html/* ./
git add --all
git commit -m "automatic commit"
git push
