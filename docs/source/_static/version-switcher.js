document.addEventListener("DOMContentLoaded", async () => {
  const marker = "/docs/";
  const markerIndex = window.location.pathname.indexOf(marker);
  const sidebar = document.querySelector(".wy-side-nav-search");
  if (markerIndex < 0 || !sidebar) return;

  const siteRoot = window.location.pathname.slice(0, markerIndex + 1);
  const currentVersion = window.location.pathname
    .slice(markerIndex + marker.length)
    .split("/")[0];

  try {
    const response = await fetch(`${siteRoot}versions.json`);
    if (!response.ok) return;
    const versions = await response.json();
    if (!versions.length) return;

    const container = document.createElement("div");
    container.className = "documentation-version";

    const label = document.createElement("label");
    label.htmlFor = "documentation-version-select";
    label.textContent = "Documentation version";

    const select = document.createElement("select");
    select.id = label.htmlFor;
    for (const version of versions) {
      const option = document.createElement("option");
      option.value = version.slug;
      option.textContent =
        version.kind === "pull_request"
          ? `${version.label}: ${version.title}`
          : version.label;
      option.selected = version.slug === currentVersion;
      select.appendChild(option);
    }
    select.addEventListener("change", () => {
      window.location.href = `${siteRoot}docs/${select.value}/`;
    });

    const links = document.createElement("div");
    links.className = "documentation-version-links";

    const allVersions = document.createElement("a");
    allVersions.href = `${siteRoot}versions.html`;
    allVersions.textContent = "All versions";
    links.appendChild(allVersions);

    const current = versions.find(version => version.slug === currentVersion);
    if (current?.coverage) {
      const coverage = document.createElement("a");
      coverage.href = `${siteRoot}docs/${currentVersion}/cov/`;
      coverage.textContent = "Coverage";
      links.appendChild(coverage);
    }

    container.append(label, select, links);
    sidebar.appendChild(container);
  } catch {
    // Local and partial documentation builds do not have versions.json.
  }
});
