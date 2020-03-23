# Convert gif animations to the webp format
find . -name "*.gif" -type f -exec gif2webp {} -o {}.webp \;
find . -name "*.webp" -exec rename 's/(.*)\.gif\.webp/$1.webp/' {} \;
find . -name "*.gif" -type f -delete
find . -name "*.gif.webp" -type f -delete  # in case of duplicates
