resource "azurerm_resource_group" "rg" {
  name     = var.rg
  location = var.region

  tags = var.tags
}

data "azurerm_resource_group" "rg" {
  name       = var.rg
  depends_on = [azurerm_resource_group.rg]
}
