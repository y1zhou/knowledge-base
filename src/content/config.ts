import { z, defineCollection } from "astro:content";

const postsCollection = defineCollection({
  type: "content", // v2.5.0 and later
  schema: z.object({
    title: z.string(),
    pubDate: z.union([z.date(), z.string().transform((str) => new Date(str))]),
    author: z.string(),
    description: z.string().optional(),
    tags: z.array(z.string()).optional(),
    cover: z
      .object({
        url: z.string(),
        square: z.string(),
        alt: z.string(),
      })
      .optional(),
    theme: z.enum(["light", "dark"]).optional(),
    featured: z.boolean().optional(),
    meta: z
      .array(
        z.object({
          name: z.string(),
          content: z.string(),
        })
      )
      .optional(),
    keywords: z.string().optional(),
  }),
});

export const collections = {
  posts: postsCollection,
};
