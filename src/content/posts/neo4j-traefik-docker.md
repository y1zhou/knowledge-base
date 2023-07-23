---
title: "Neo4j Bolt behind traefik in Docker"

pubDate: 2022-05-08
description: ""
author: "Yi Zhou"
cover:
  url: "https://pic.lookcos.cn/i/usr/uploads/2022/04/2067928922.png"
  square: "https://pic.lookcos.cn/i/usr/uploads/2022/04/2067928922.png"
  alt: "cover"
tags: ["DevOps", "Web"]
# theme: 'light'
# featured: true
keywords: DevOps, Neo4j, Traefik, Docker
---

## Background

[Neo4j](https://neo4j.com/) is a graph database that provides performant queries and data science analytic pipelines. The [official manual](https://neo4j.com/docs/operations-manual/current/docker/) on running Neo4j with Docker covers all Neo4j-related configurations, but doesn’t talk about the process of putting Neo4j behind a reverse proxy.

[Traefik](https://traefik.io/) is a simple reverse proxy and load balancer that can automatically discover services within Docker. I was a Nginx user for years, but editing the config file for every container gets tiring pretty quickly. With traefik, we just need to define some extra labels while writing the [Docker Compose](https://docs.docker.com/compose/) file.

## Prerequisites

- A running [Docker](https://docs.docker.com/get-started/overview/) daemon.
- The [Compose](https://github.com/docker/compose) plugin for Docker.
- Access to a domain name.

## Compose files

### Traefik

First we need to get traefik up and running. Note that my domain is registered with CloudFlare, so the SSL configuration part uses `cloudflare` as its provider. You may of course use other providers as suggested by [traefik’s documentation](https://doc.traefik.io/traefik/v1.6/configuration/acme/#provider).

In your terminal, `cd` into a directory for storing traefik-related files (I like `~/docker`), and run the following commands:

```bash
mkdir -p ~/docker/traefik && cd ~/docker/traefik
mkdir logs
mkdir letsencrypt
touch docker-compose.yaml
```

These create two directories named `logs/` and `letsencrypt/`, and an empty file called `docker-compose.yaml`. Now use your favorite text editor to open `docker-compose-yaml` and paste in the following:

```yaml
version: "3"

services:
  traefik:
    image: "traefik:v2.6"
    container_name: "traefik"
    restart: always
    networks:
      - traefik
    command:
      # Docker configuration
      - "--providers.docker=true"
      - "--providers.docker.exposedbydefault=false"

      # Logging
      - "--accesslog.filepath=/logs/access.log"
      - "--accesslog.bufferingsize=100"

      # Entrypoint
      - "--entrypoints.web.address=:80"
      - "--entrypoints.websecure.address=:443"

      # SSL configuration
      - "--certificatesresolvers.myresolver.acme.dnschallenge=true"
      - "--certificatesresolvers.myresolver.acme.dnschallenge.provider=cloudflare"
      - "--certificatesresolvers.myresolver.acme.email=me@example.com"
      - "--certificatesresolvers.myresolver.acme.storage=/letsencrypt/acme.json"
      # HTTP -> HTTPS
      - "--entrypoints.web.http.redirections.entryPoint.to=websecure"
      - "--entrypoints.web.http.redirections.entryPoint.scheme=https"

    environment:
      - "CF_DNS_API_TOKEN=@@@my-secret-token@@@"
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - "./letsencrypt:/letsencrypt"
      - "/var/run/docker.sock:/var/run/docker.sock:ro"
      - "./logs:/logs"

networks:
  traefik:
    name: traefik
```

Make sure you have the correct API token and email filled in. Now to spin up traefik, simply run:

```bash
docker compose up -d
```

### Neo4j

Next we have the interesting part. The Neo4j Docker images comes with the database server and a built-in Neo4j Browser. The browser (and other API drivers) communicates with the server using the Bolt protocol. This means to have a fully functional Neo4j behind traefik, the proxy needs to handle HTTP(S) traffic to the browser endpoint _and_ TCP traffic in Bolt.

<aside>
🔥 You don’t really need to worry about the Bolt part if you are only going to use the Neo4j Browser to interact with the database. The browser talks to the backend through a WebSocket.

</aside>

After `cd`-ing into the working directory for neo4j, we first create an `.env` file containing the default username and password for the admin user.

```bash
mkdir -p ~/docker/neo4j && cd ~/docker/neo4j
echo 'NEO4J_AUTH=neo4j/my-password' > neo4j.env
```

Change `my-password` to something more secure. Next we use `id -u` and `id -g` to find the UID and GID of the user. These values will be used to replace `services.neo4j.user` in the YAML file below. Here’s the `docker-compose` file for the Neo4j container:

```yaml
version: "3"

services:
  neo4j:
    image: "neo4j:4.4.6" # or "neo4j:4.4.6-enterprise"
    restart: unless-stopped
    user: "1042:1042" # use `id -u` and `id -g` to find UID and GID
    ports:
      - "7687:7687"
    volumes:
      - neo4j-data:/data
    env_file:
      - neo4j.env
    environment:
      - NEO4J_dbms_default__listen__address=0.0.0.0
      - NEO4J_dbms_connector_bolt_advertised__address=n4j.example.com:443
      - NEO4J_dbms_allow__upgrade=true
      - NEO4J_ACCEPT_LICENSE_AGREEMENT=yes # if using enterprise edition
      - NEO4J_metrics_enabled=false # a little bit of privacy

    networks:
      - traefik
    labels:
      - "traefik.enable=true"

      # Neo4j browser service
      - "traefik.http.routers.neo4j.rule=Host(`n4j.example.com`) && PathPrefix(`/neo4j`) || PathPrefix(`/browser`)"
      - "traefik.http.routers.neo4j.service=neo4j"
      - "traefik.http.services.neo4j.loadbalancer.server.port=7474"
      - "traefik.http.routers.neo4j.entrypoints=websecure"
      - "traefik.http.routers.neo4j.tls.certresolver=myresolver"

      # neo4j bolt service for browser websocket
      - "traefik.http.routers.neo4j-bolt.rule=Host(`n4j.example.com`)"
      - "traefik.http.routers.neo4j-bolt.service=neo4j-bolt"
      - "traefik.http.routers.neo4j-bolt.entrypoints=websecure"
      - "traefik.http.services.neo4j-bolt.loadbalancer.server.port=7687"
      - "traefik.http.routers.neo4j-bolt.tls.certresolver=myresolver"

      - "traefik.http.middlewares.sslheader.headers.customrequestheaders.X-Forwarded-Proto=https,wss"
      - "traefik.http.routers.neo4j-bolt.middlewares=sslheader"

      # Bolt service for drivers
      - "traefik.tcp.routers.neo4j-bolt.rule=HostSNI(`n4j.example.com`)"
      - "traefik.tcp.routers.neo4j-bolt.service=neo4j-bolt"
      - "traefik.tcp.routers.neo4j-bolt.entrypoints=websecure"
      - "traefik.tcp.routers.neo4j-bolt.tls.passthrough=true"
      - "traefik.tcp.services.neo4j-bolt.loadbalancer.server.port=7687"
      - "traefik.tcp.routers.neo4j-bolt.tls.certresolver=myresolver"

volumes:
  neo4j-data:

networks:
  traefik:
    external: true
    name: traefik
```

Replace `n4j.example.com` with your desired subdomain. Now run `docker compose up -d`, wait a minute, and open `n4j.example.com/browser` to get started.

![Neo4j Browser interface for connecting to the database.](~/assets/images/neo4j-browser-interface.png)

Neo4j Browser interface for connecting to the database.

Note that in “Connect URL” we are using port 443 instead of port 7687. Type in your predefined username and password from the `.env` file, hit connect, and you can start messing around with the graph!

As for connecting to the database with [drivers](https://neo4j.com/developer/language-guides/), we use port 7687 and drop the `+s` in the schema. Taking the Python driver as an example:

```python
from neo4j import GraphDatabase
conn = GraphDatabase.driver(
  'neo4j://n4j.example.com',  # :7687 is used by default
  auth=('neo4j_username', 'neo4j_password')
)
conn.verify_connectivity()
```

Not the safest implementation, but I haven’t figured out a way to get Neo4j to use the auto-generated certificates from traefik yet. Note that this is connecting to your server’s 7687 port directly, so **make sure you are using a strong password**!
