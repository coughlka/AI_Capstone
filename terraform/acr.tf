# Azure Container Registry
resource "azurerm_container_registry" "acr" {
  name                = var.acr_name
  resource_group_name = azurerm_resource_group.rg.name
  location            = var.region
  sku                 = "Basic"
  admin_enabled       = true

  tags = var.tags
}

# =============================================================================
# ACR Webhook — triggers App Service to pull the latest image after CI/CD push
# =============================================================================
resource "azurerm_container_registry_webhook" "webapp_deploy" {
  name                = "webappDeploy"
  registry_name       = azurerm_container_registry.acr.name
  resource_group_name = azurerm_resource_group.rg.name
  location            = var.region

  service_uri = "https://${azurerm_linux_web_app.app_service.site_credential[0].name}:${azurerm_linux_web_app.app_service.site_credential[0].password}@${azurerm_linux_web_app.app_service.name}.scm.azurewebsites.net/api/registry/webhook"

  actions = ["push"]
  status  = "enabled"
  scope   = "biomarker-cancer-web-app:latest"

  tags = var.tags
}

# =============================================================================
# Outputs — use these to set GitHub secrets
# =============================================================================
output "acr_login_server" {
  description = "ACR login server URL"
  value       = azurerm_container_registry.acr.login_server
}

output "acr_admin_username" {
  description = "ACR admin username"
  value       = azurerm_container_registry.acr.admin_username
}

output "acr_admin_password" {
  description = "ACR admin password"
  value       = azurerm_container_registry.acr.admin_password
  sensitive   = true
}
