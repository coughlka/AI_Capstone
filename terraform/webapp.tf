# Azure Web App
resource "azurerm_service_plan" "app_service_plan" {
  name                = var.app_service_plan_name
  resource_group_name = azurerm_resource_group.rg.name
  location            = var.region
  os_type             = "Linux"
  sku_name            = var.sku_size

  tags = var.tags
}

resource "azurerm_linux_web_app" "app_service" {
  name                = var.app_service_name
  resource_group_name = azurerm_resource_group.rg.name
  location            = var.region
  service_plan_id     = azurerm_service_plan.app_service_plan.id

  site_config {
    always_on              = true
    vnet_route_all_enabled = true
  }

  virtual_network_subnet_id = azurerm_subnet.webapp_subnet.id

  storage_account {
    name         = var.storage_account_name
    share_name   = azurerm_storage_share.web_app_share.name
    type         = "AzureFiles"
    access_key   = azurerm_storage_account.storage_account.primary_access_key
    account_name = azurerm_storage_account.storage_account.name
    mount_path   = "/web-app"
  }

  tags = var.tags
}
