#!/usr/bin/env bash

# Convert gif animations to the webp format
find ./content -name "*.gif" -type f -exec gif2webp {} -o {}.webp \;
find ./content -name "*.webp" -exec rename 's/(.*)\.gif\.webp$/$1.webp/' {} \;
find ./content -name "*.gif" -type f -delete
find ./content -name "*.gif.webp" -type f -delete  # in case of duplicates
