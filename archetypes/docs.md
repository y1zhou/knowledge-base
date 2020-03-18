---
title: "{{ replace .Name "-" " " | title }}"
slug: "{{ .Name | title | lower }}"
categories:
  - Statistics
tags:
  - Statistics

summary: 
date: {{ .Date }}
toc: true
type: docs  # Do not modify.
weight: 1000

menu:
  example:
    name: YourPageTitle
    # identifier: YourPageID
    # parent: YourParentID
    weight: 1
---
