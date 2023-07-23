## 👏🏻 Introduction

This is a minimalist, beautiful, responsive blogging program written in Astro.

### Three display model of images

The three display modes of images are: `wide`, `big`, `inline`.
When you edit your markdown file, you can add `wide` or `big` or `inline` to the image alt, like this:

```markdown
![alt content|wide](a.png)
```

<strong>The Separator is `|`, and the default mode is `big`.</strong>

## 🚀 Project Structure

Astro looks for `.astro` or `.md` files in the `src/content/` directory. Each page is exposed as a route based on its file name. This is known as [content collections](https://docs.astro.build/en/guides/content-collections/).

There's nothing special about `src/components/`, but that's where we like to put any Astro/React/Vue/Svelte/Preact [components](https://docs.astro.build/en/core-concepts/framework-components/).

Any static assets, like [images](https://docs.astro.build/en/guides/images/), can be placed in the `src/assets/images` directory.

## 🧞 Commands

All commands are run from the root of the project, from a terminal:

| Command                | Action                                           |
| :--------------------- | :----------------------------------------------- |
| `npm install`          | Installs dependencies                            |
| `npm run dev`          | Starts local dev server at `localhost:3000`      |
| `npm run build`        | Build your production site to `./dist/`          |
| `npm run preview`      | Preview your build locally, before deploying     |
| `npm run astro ...`    | Run CLI commands like `astro add`, `astro check` |
| `npm run astro --help` | Get help using the Astro CLI                     |

## 👀 Preview

[https://yufengbiji.com](https://yufengbiji.com)

[https://astro.yufengbiji.com](https://astro.yufengbiji.com)
